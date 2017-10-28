### List of Ensembl MySQL servers / ports
###   https://www.ensembl.org/info/data/mysql.html
### Ensembl Core schema
###   https://www.ensembl.org/info/docs/api/core/core_schema.html

library(GenomicFeatures)
library(RMySQL)

.fix_numeric_cols <- function(df)
{
    col_idx <- which(sapply(df, is.numeric))
    df[col_idx] <- lapply(df[col_idx], as.integer)
    df
}

.RMySQL_SELECT <- function(dbconn, SQL)
{
    ## Not sure systematic conversion of numeric to int is actually a
    ## good idea (risk of overflow?)
    .fix_numeric_cols(suppressWarnings(dbGetQuery(dbconn, SQL)))
}

.seq_region_columns <- c(
    "seq_region_id",
    "seq_region_start",
    "seq_region_end",
    "seq_region_strand"
)

.get_transcripts <- function(dbconn)
{
    transcript_columns <- c(
        "transcript_id",
        "stable_id",
        .seq_region_columns,
        "biotype"
    )
    SELECT_columns <- c(paste0("transcript.", transcript_columns),
                        "gene.stable_id")
    SQL <- sprintf("SELECT %s FROM transcript LEFT JOIN gene USING(gene_id)",
                   paste(SELECT_columns, collapse=","))
    transcripts <- .RMySQL_SELECT(dbconn, SQL)
    colnames(transcripts) <- c("tx_id",
                               "tx_name",
                               "seq_region_id",
                               "tx_start",
                               "tx_end",
                               "tx_strand",
                               "tx_type",
                               "gene_id")
    transcripts$tx_strand <- strand(transcripts$tx_strand)
    transcripts
}

.get_splicings <- function(dbconn)
{
    exon_columns <- c(
        "exon_id",
        "stable_id",
        .seq_region_columns
    )
    SELECT_columns <- c("transcript_id", "rank", exon_columns)
    SQL <- sprintf("SELECT %s FROM exon_transcript INNER JOIN exon USING(exon_id)",
                   paste(SELECT_columns, collapse=","))
    splicings <- .RMySQL_SELECT(dbconn, SQL)
    colnames(splicings) <- c("tx_id",
                             "exon_rank",
                             "exon_id",
                             "exon_name",
                             "seq_region_id",
                             "exon_start",
                             "exon_end",
                             "exon_strand")
    splicings$exon_strand <- strand(splicings$exon_strand)
    splicings
}

.get_cds <- function(dbconn)
{
    translation_columns <- c(
        "translation_id",
        "stable_id",
        "start_exon_id",   # ==> exon.exon_id
        "seq_start",       # relative to first exon
        "end_exon_id",     # ==> exon.exon_id
        "seq_end",         # relative to last exon
        "transcript_id"
    )
    SQL <- sprintf("SELECT %s FROM translation",
                   paste(translation_columns, collapse=","))
    .RMySQL_SELECT(dbconn, SQL)
}

.get_toplevel_seq_region_ids <- function(dbconn)
{
    ## FIXME: id0 is hard-coded to 6 for now.
    ## See .Ensembl_fetchAttribTypeIdForTopLevelSequence() for how to
    ## extract this value from Ensembl
    id0 <- 6L
    seq_region_attrib_columns <- c(
        "seq_region_id",
        "attrib_type_id",
        "value"
    )
    SQL <- sprintf("SELECT %s FROM seq_region_attrib",
                   paste(seq_region_attrib_columns, collapse=","))
    seq_region_attrib <- .RMySQL_SELECT(dbconn, SQL)
    seq_region_attrib$seq_region_id[seq_region_attrib$attrib_type_id == id0]
}

.get_chromlengths <- function(dbconn, seq_region_ids=NULL)
{
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
    SQL <- "SELECT * FROM seq_region INNER JOIN coord_system USING(coord_system_id)"
    seq_region <- .RMySQL_SELECT(dbconn, SQL)
    stopifnot(identical(colnames(seq_region), joined_columns))
    colnames(seq_region)[6:9] <- paste0("coord_system_",
                                        colnames(seq_region)[6:9])
    top_level_ids <- .get_toplevel_seq_region_ids(dbconn)
    chromlengths <- GenomicFeatures:::extract_chromlengths_from_seq_region(
                                              seq_region,
                                              top_level_ids,
                                              seq_region_ids=seq_region_ids)
    colnames(chromlengths)[2L] <- "chrom"
    chromlengths
}

makeTxDbFromEnsembl <- function(dbname="homo_sapiens_core_90_38",
                                host="useastdb.ensembl.org")
{
    dbconn <- dbConnect(MySQL(), dbname=dbname, username="anonymous", host=host)
    on.exit(dbDisconnect(dbconn))

    transcripts <- .get_transcripts(dbconn)

    genes <- transcripts[ , c("tx_name", "gene_id")]
    transcripts$gene_id <- NULL

    splicings <- .get_splicings(dbconn)

    #cds <- .get_cds(dbconn)

    seq_region_ids <- unique(c(transcripts$seq_region_id,
                               splicings$seq_region_id))

    chromlengths <- .get_chromlengths(dbconn, seq_region_ids=seq_region_ids)

    m <- match(transcripts$seq_region_id, chromlengths$seq_region_id)
    transcripts$tx_chrom <- chromlengths$chrom[m]
    transcripts$seq_region_id <- NULL

    m <- match(splicings$seq_region_id, chromlengths$seq_region_id)
    splicings$exon_chrom <- chromlengths$chrom[m]
    splicings$seq_region_id <- NULL

    chromlengths$seq_region_id <- NULL

    makeTxDb(transcripts, splicings, genes=genes, chrominfo=chromlengths,
             reassign.ids=TRUE)
}

