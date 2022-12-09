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
    mysql_url <- ftp_url_to_Ensembl_mysql(release)
    core_dirs <- Ensembl_listMySQLCoreDirs(mysql_url, release=release)
    prefix <- paste0(gsub(" ", "_", tolower(organism), fixed=TRUE), "_core_")
    i <- match(prefix, substr(core_dirs, 1L, nchar(prefix)))
    if (is.na(i)) {
        pattern <- "^([^_])[^_]*_(.*_core_.*)$"
        replacement <- "\\1\\2"
        alt_core_dirs <- sub(pattern, replacement, core_dirs)
        i <- match(prefix, substr(alt_core_dirs, 1L, nchar(prefix)))
        if (is.na(i)) {
            if (is.na(release))
                where <- "the current Ensembl release"
            else
                where <- paste0("Ensembl release ", release)
            stop(wmsg("no core db found in ", where, " ",
                      "for organism \"", organism, "\""))
        }
    }
    dbname <- core_dirs[[i]]
    if (!is.na(release))  # sanity check
        stopifnot(.dbname2release(dbname) == release)
    dbname
}

.dbselect <- function(dbconn, columns, from)
{
    SQL <- sprintf("SELECT %s FROM %s", paste0(columns, collapse=","), from)
    dbGetQuery(dbconn, SQL)
}

.fetch_attrib_type_id <- function(dbconn, tx_attrib)
{
    from <- sprintf("attrib_type WHERE code='%s'", tx_attrib)
    attrib_type_id <- .dbselect(dbconn, "attrib_type_id", from)[[1L]]
    if (length(attrib_type_id) == 0L)
        stop(wmsg("invalid supplied transcript attribute code: ", tx_attrib))
    attrib_type_id
}

.restrict_to_tx_attrib <- function(from, tx_attrib_type_id)
{
    paste0(from, " JOIN transcript_attrib AS ta",
                 " USING(transcript_id)",
                 " WHERE ta.attrib_type_id='", tx_attrib_type_id, "'")
}

.seq_region_columns <- c(
    "seq_region_id",
    "seq_region_start",
    "seq_region_end",
    "seq_region_strand"
)

.fetch_Ensembl_transcripts <- function(dbconn, tx_attrib_type_id=NULL)
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
    if (!is.null(tx_attrib_type_id))
        from <- .restrict_to_tx_attrib(from, tx_attrib_type_id)
    transcripts <- .dbselect(dbconn, columns, from)
    colnames(transcripts) <- c("tx_id",
                               "tx_name",
                               "seq_region_id",
                               "tx_start",
                               "tx_end",
                               "tx_strand",
                               "tx_type",
                               "gene_id")
    transcripts$tx_strand <- strand(transcripts$tx_strand)
    nb_transcripts <- length(unique(transcripts$tx_id))
    nb_genes <- length(unique(transcripts$gene_id))
    message("OK")
    message("  (fetched ", nb_transcripts, " transcripts from ",
            nb_genes, " genes)")
    transcripts
}

.fetch_Ensembl_translations <- function(dbconn, tx_attrib_type_id=NULL)
{
    columns <- c(
        "stable_id",
        "start_exon_id",   # ==> exon.exon_id
        "seq_start",       # relative to first exon
        "end_exon_id",     # ==> exon.exon_id
        "seq_end",         # relative to last exon
        "transcript_id"
    )
    from <- "translation"
    if (!is.null(tx_attrib_type_id))
        from <- .restrict_to_tx_attrib(from, tx_attrib_type_id)
    .dbselect(dbconn, columns, from)
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
.add_cds_cols <- function(dbconn, splicings, tx_attrib_type_id=NULL)
{
    translations <- .fetch_Ensembl_translations(dbconn, tx_attrib_type_id)

    m <- match(splicings$tx_id, translations$transcript_id)
    cds_name <- translations$stable_id[m]

    m1 <- matchIntegerPairs(splicings$tx_id,
                                        splicings$exon_id,
                                        translations$transcript_id,
                                        translations$start_exon_id)
    offset1 <- translations$seq_start[m1] - 1L
    m2 <- matchIntegerPairs(splicings$tx_id,
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

    cbind(splicings, cds_name=cds_name, cds_start=cds_start, cds_end=cds_end,
                     stringsAsFactors=FALSE)
}

.fetch_Ensembl_splicings <- function(dbconn, tx_attrib_type_id=NULL)
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
    if (!is.null(tx_attrib_type_id))
        from <- .restrict_to_tx_attrib(from, tx_attrib_type_id)
    splicings <- .dbselect(dbconn, columns, from)
    colnames(splicings) <- c("tx_id",
                             "exon_rank",
                             "exon_id",
                             "exon_name",
                             "seq_region_id",
                             "exon_start",
                             "exon_end",
                             "exon_strand")
    splicings$exon_strand <- strand(splicings$exon_strand)
    oo <- orderIntegerPairs(splicings$tx_id, splicings$exon_rank)
    splicings <- S4Vectors:::extract_data_frame_rows(splicings, oo)
    splicings <- .add_cds_cols(dbconn, splicings, tx_attrib_type_id)
    message("OK")
    splicings
}

.get_toplevel_seq_region_ids <- function(dbconn)
{
    ## FIXME: id0 is hard-coded to 6 for now.
    ## See .Ensembl_fetchAttribTypeIdForTopLevelSequence() for how
    ## to extract this value from Ensembl.
    id0 <- 6L
    columns <- c("seq_region_id", "attrib_type_id", "value")
    seq_region_attrib <- .dbselect(dbconn, columns, "seq_region_attrib")
    seq_region_attrib$seq_region_id[seq_region_attrib$attrib_type_id == id0]
}

.fetch_Ensembl_chrominfo <- function(dbconn, seq_region_ids=NULL,
                                     circ_seqs=NULL)
{
    message("Fetch chromosome names and lengths from Ensembl ...",
            appendLF=FALSE)
    seq_region_columns <- c(
        "seq_region_id",
        "seq_region.name AS name",
        "coord_system_id",
        "length"
    )
    coord_system_columns <- c(
        "coord_system_id",
        "species_id",
        "coord_system.name AS coord_system_name",
        "coord_system.version AS coord_system_version",
        "coord_system.rank AS coord_system_rank",
        "coord_system.attrib AS coord_system_attrib"
    )
    using_column <- "coord_system_id"
    joined_columns <- c(using_column,
                        setdiff(seq_region_columns, using_column),
                        setdiff(coord_system_columns, using_column))
    from <- "seq_region INNER JOIN coord_system USING(coord_system_id)"
    seq_region <- .dbselect(dbconn, joined_columns, from)
    top_level_ids <- .get_toplevel_seq_region_ids(dbconn)
    chromlengths <- extract_chromlengths_from_seq_region(
                                         seq_region,
                                         top_level_ids,
                                         seq_region_ids=seq_region_ids)
    chrominfo <- data.frame(
        seq_region_id=chromlengths$seq_region_id,
        chrom=chromlengths$name,
        length=chromlengths$length,
        is_circular=GenomeInfoDb:::make_circ_flags_from_circ_seqs(
                                                   chromlengths$name,
                                                   circ_seqs)
    )
    message("OK")
    chrominfo
}

.gather_Ensembl_metadata <- function(organism, dbname, server,
                                     tx_attrib=NULL)
{
    message("Gather the metadata ... ", appendLF=FALSE)
    release <- .dbname2release(dbname)
    full_dataset <- is.null(tx_attrib)
    metadata <- data.frame(name=c("Data source",
                                  "Organism",
                                  "Ensembl release",
                                  "Ensembl database",
                                  "MySQL server",
                                  "Full dataset"),
                           value=c("Ensembl",
                                   organism,
                                   release,
                                   dbname,
                                   server,
                                   ifelse(full_dataset, "yes", "no")),
                           stringsAsFactors=FALSE)
    if (!full_dataset) {
        val <- sprintf("transcripts with attribute \"%s\"", tx_attrib)
        more_metadata <- data.frame(name="Imported only", value=val,
                                    stringsAsFactors=FALSE)
        metadata <- rbind(metadata, more_metadata)
    }
    message("OK")
    metadata
}

### Always set 'server' to "useastdb.ensembl.org" in the examples so that
### they run fast on the build machines (which are located on the East Coast).
makeTxDbFromEnsembl <- function(organism="Homo sapiens",
                                release=NA,
                                circ_seqs=NULL,
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

    if (is.null(tx_attrib)) {
        tx_attrib_type_id <- NULL
    } else {
        if (!isSingleString(tx_attrib))
            stop(wmsg("'tx_attrib' must be a single string or NULL"))
        tx_attrib_type_id <- .fetch_attrib_type_id(dbconn, tx_attrib)
    }

    transcripts <- .fetch_Ensembl_transcripts(dbconn, tx_attrib_type_id)

    splicings <- .fetch_Ensembl_splicings(dbconn, tx_attrib_type_id)

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

    metadata <- .gather_Ensembl_metadata(organism, dbname, server, tx_attrib)

    message("Make the TxDb object ... ", appendLF=FALSE)
    txdb <- makeTxDb(transcripts, splicings,
                     chrominfo=chrominfo, metadata=metadata, reassign.ids=TRUE)
    message("OK")
    txdb
}

