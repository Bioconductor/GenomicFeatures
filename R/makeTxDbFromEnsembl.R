### =========================================================================
### makeTxDbFromEnsembl()
### -------------------------------------------------------------------------
###
### List of Ensembl public MySQL servers / ports
###   https://www.ensembl.org/info/data/mysql.html
### Ensembl Core schema
###   https://www.ensembl.org/info/docs/api/core/core_schema.html
###


.normarg_organism <- function(organism)
{
    if (!isSingleString(organism))
        stop("'organism' must be a single string")
    ## Remove extra spaces.
    tmp <- strsplit(organism, " ")[[1L]]
    paste(tmp[nzchar(tmp)], collapse=" ")
}

.dbname2release <- function(dbname)
    as.integer(gsub("^.*_core_([0-9]+)_.*$", "\\1", dbname))

.lookup_dbname <- function(organism, release=NA)
{
    organism <- .normarg_organism(organism)
    if (!isSingleNumberOrNA(release))
        stop("'release' must be a valid Ensembl release number or NA")
    available_dbs <- Ensembl_listMySQLCoreDirs(release=release)
    prefix <- paste0(gsub(" ", "_", tolower(organism), fixed=TRUE), "_core_")
    i <- match(prefix, substr(available_dbs, 1L, nchar(prefix)))
    dbname <- available_dbs[[i]]
    if (!is.na(release))  # sanity check
        stopifnot(.dbname2release(dbname) == release)
    dbname
}

.restrict_to_tx_attrib <- function(SQL, tx_attrib) {
    paste0(SQL,
           " JOIN transcript_attrib AS ta",
           " USING(transcript_id)",
           " WHERE ta.attrib_type_id=",
           "(SELECT attrib_type_id FROM attrib_type WHERE code='", tx_attrib,
           "')")
}

.dbselect <- function(dbconn, columns, from, tx_attrib = NULL)
{
    SQL <- sprintf("SELECT %s FROM %s", paste0(columns, collapse=","),
                   from)
    if (!is.null(tx_attrib)) {
        SQL <- .restrict_to_tx_attrib(SQL, tx_attrib)
    }
    dbGetQuery(dbconn, SQL)
}

.seq_region_columns <- c(
    "seq_region_id",
    "seq_region_start",
    "seq_region_end",
    "seq_region_strand"
)

.fetch_Ensembl_transcripts <- function(dbconn, tx_attrib)
{
    message("Fetch transcripts and genes from Ensembl ... ",
            appendLF=FALSE)
    transcript_columns <- c(
        "transcript_id",
        "stable_id",
        .seq_region_columns,
        "biotype"
    )
    columns <- c(paste0("transcript.", transcript_columns), "gene.stable_id")
    from <- "transcript LEFT JOIN gene USING(gene_id)"
    transcripts <- .dbselect(dbconn, columns, from, tx_attrib)
    colnames(transcripts) <- c("tx_id",
                               "tx_name",
                               "seq_region_id",
                               "tx_start",
                               "tx_end",
                               "tx_strand",
                               "tx_type",
                               "gene_id")
    transcripts$tx_strand <- strand(transcripts$tx_strand)
    message("OK")
    transcripts
}

.fetch_Ensembl_translations <- function(dbconn, tx_attrib)
{
    columns <- c(
        "stable_id",
        "start_exon_id",   # ==> exon.exon_id
        "seq_start",       # relative to first exon
        "end_exon_id",     # ==> exon.exon_id
        "seq_end",         # relative to last exon
        "transcript_id"
    )
    .dbselect(dbconn, columns, "translation", tx_attrib)
}

### 'has_cds' must be a logical vector.
### 'tx_id' must be an atomic vector parallel to 'has_cds'.
.has_cds <- function(has_cds, tx_id)
{
    stopifnot(length(has_cds) == length(tx_id))
    breakpoints <- cumsum(runLength(Rle(tx_id)))
    partitioning <- PartitioningByEnd(breakpoints)
    ## List of relative CDS indices.
    ridx_list <- which(relist(has_cds, partitioning))
    idx0 <- which(lengths(ridx_list) == 0L)
    min_ridx <- replaceROWS(min(ridx_list), idx0, 1L)
    max_ridx <- replaceROWS(max(ridx_list), idx0, 0L)
    ## Absolute CDS index.
    aidx <- shift(IRanges(min_ridx, max_ridx), start(partitioning) - 1L)
    has_cds <- logical(length(tx_id))
    replaceROWS(has_cds, aidx, TRUE)
}

### Add "cds_name", "cds_start", and "cds_end" cols to 'splicings'.
.add_cds_cols <- function(dbconn, splicings, tx_attrib)
{
    translations <- .fetch_Ensembl_translations(dbconn, tx_attrib)

    m <- match(splicings$tx_id, translations$transcript_id)
    cds_name <- translations$stable_id[m]

    m1 <- S4Vectors:::matchIntegerPairs(splicings$tx_id,
                                        splicings$exon_id,
                                        translations$transcript_id,
                                        translations$start_exon_id)
    offset1 <- translations$seq_start[m1] - 1L
    m2 <- S4Vectors:::matchIntegerPairs(splicings$tx_id,
                                        splicings$exon_id,
                                        translations$transcript_id,
                                        translations$end_exon_id)
    offset2 <- translations$seq_end[m2] - 1L

    cds_start <- ifelse(splicings$exon_strand == "+",
                        splicings$exon_start + offset1,
                        splicings$exon_end - offset2)
    cds_end <- ifelse(splicings$exon_strand == "+",
                      splicings$exon_start + offset2,
                      splicings$exon_end - offset1)

    has_cds <- .has_cds(!(is.na(cds_start) & is.na(cds_end)), splicings$tx_id)
    cds_name[!has_cds] <- NA
    idx <- which(has_cds & is.na(cds_start))
    cds_start[idx] <- splicings$exon_start[idx]
    idx <- which(has_cds & is.na(cds_end))
    cds_end[idx] <- splicings$exon_end[idx]

    cbind(splicings, cds_name, cds_start, cds_end, stringsAsFactors=FALSE)
}

.fetch_Ensembl_splicings <- function(dbconn, tx_attrib)
{
    message("Fetch exons and CDS from Ensembl ... ",
            appendLF=FALSE)
    exon_columns <- c(
        "exon_id",
        "stable_id",
        .seq_region_columns
    )
    columns <- c("transcript_id", "rank", exon_columns)
    from <- "exon_transcript INNER JOIN exon USING(exon_id)"
    splicings <- .dbselect(dbconn, columns, from, tx_attrib)
    colnames(splicings) <- c("tx_id",
                             "exon_rank",
                             "exon_id",
                             "exon_name",
                             "seq_region_id",
                             "exon_start",
                             "exon_end",
                             "exon_strand")
    splicings$exon_strand <- strand(splicings$exon_strand)
    oo <- S4Vectors:::orderIntegerPairs(splicings$tx_id, splicings$exon_rank)
    splicings <- S4Vectors:::extract_data_frame_rows(splicings, oo)
    splicings <- .add_cds_cols(dbconn, splicings, tx_attrib)
    message("OK")
    splicings
}

.get_toplevel_seq_region_ids <- function(dbconn)
{
    ## FIXME: id0 is hard-coded to 6 for now.
    ## See .Ensembl_fetchAttribTypeIdForTopLevelSequence() for how to
    ## extract this value from Ensembl
    id0 <- 6L
    columns <- c("seq_region_id", "attrib_type_id", "value")
    seq_region_attrib <- .dbselect(dbconn, columns, "seq_region_attrib")
    seq_region_attrib$seq_region_id[seq_region_attrib$attrib_type_id == id0]
}

.fetch_Ensembl_chrominfo <- function(dbconn, seq_region_ids=NULL,
                                     circ_seqs=DEFAULT_CIRC_SEQS)
{
    message("Fetch chromosome names and lengths from Ensembl ...",
            appendLF=FALSE)
    seq_region_columns <- c(
        "seq_region_id",
        "name",
        "coord_system_id",
        "length"
    )
    coord_system_columns <- c(
        "coord_system_id",
        "species_id",
        "name",
        "version",
        "rank",
        "attrib"
    )
    using_column <- "coord_system_id"
    joined_columns <- c(using_column,
                        setdiff(seq_region_columns, using_column),
                        setdiff(coord_system_columns, using_column))
    from <- "seq_region INNER JOIN coord_system USING(coord_system_id)"
    seq_region <- .dbselect(dbconn, "*", from)
    stopifnot(identical(colnames(seq_region), joined_columns))
    colnames(seq_region)[6:9] <- paste0("coord_system_",
                                        colnames(seq_region)[6:9])
    top_level_ids <- .get_toplevel_seq_region_ids(dbconn)
    chromlengths <- extract_chromlengths_from_seq_region(
                                         seq_region,
                                         top_level_ids,
                                         seq_region_ids=seq_region_ids)
    chrominfo <- data.frame(
        seq_region_id=chromlengths$seq_region_id,
        chrom=chromlengths$name,
        length=chromlengths$length,
        is_circular=make_circ_flags_from_circ_seqs(chromlengths$name,
                                                   circ_seqs)
    )
    message("OK")
    chrominfo
}

.gather_Ensembl_metadata <- function(organism, dbname, server)
{
    message("Gather the metadata ... ", appendLF=FALSE)
    release <- .dbname2release(dbname)
    metadata <- data.frame(name=c("Data source",
                                  "Organism",
                                  "Ensembl release",
                                  "Ensembl database",
                                  "MySQL server"),
                           value=c("Ensembl",
                                   organism,
                                   release,
                                   dbname,
                                   server),
                           stringsAsFactors=FALSE)
    message("OK")
    metadata
}

### Always set 'server' to "useastdb.ensembl.org" in the examples so that
### they run fast on the build machines (which are located on the East Coast).
makeTxDbFromEnsembl <- function(organism="Homo sapiens",
                                release=NA,
                                circ_seqs=DEFAULT_CIRC_SEQS,
                                server="ensembldb.ensembl.org",
                                username="anonymous", password=NULL,
                                port=0L, tx_attrib=NULL)
{
    if (!requireNamespace("RMariaDB", quietly=TRUE))
        stop(wmsg("Couldn't load the RMariaDB package. ",
                  "You need to install the RMariaDB package ",
                  "in order to use makeTxDbFromEnsembl()."))

    dbname <- .lookup_dbname(organism, release=release)
    dbconn <- dbConnect(RMariaDB::MariaDB(), dbname=dbname,
                        host=server, username=username, password=password,
                        port=port)
    on.exit(dbDisconnect(dbconn))

    transcripts <- .fetch_Ensembl_transcripts(dbconn, tx_attrib)

    splicings <- .fetch_Ensembl_splicings(dbconn, tx_attrib)

    seq_region_ids <- unique(c(transcripts$seq_region_id,
                               splicings$seq_region_id))

    chrominfo <- .fetch_Ensembl_chrominfo(dbconn,
                                          seq_region_ids=seq_region_ids,
                                          circ_seqs=circ_seqs)

    m <- match(transcripts$seq_region_id, chrominfo$seq_region_id)
    transcripts$tx_chrom <- chrominfo$chrom[m]
    transcripts$seq_region_id <- NULL

    m <- match(splicings$seq_region_id, chrominfo$seq_region_id)
    splicings$exon_chrom <- chrominfo$chrom[m]
    splicings$seq_region_id <- NULL

    chrominfo$seq_region_id <- NULL

    metadata <- .gather_Ensembl_metadata(organism, dbname, server)

    message("Make the TxDb object ... ", appendLF=FALSE)
    txdb <- makeTxDb(transcripts, splicings,
                     chrominfo=chrominfo, metadata=metadata, reassign.ids=TRUE)
    message("OK")
    txdb
}

