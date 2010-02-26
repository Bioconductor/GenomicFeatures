.transcriptsORexonsORcdsBy <-
function(txdb, groupBy = c("gene", "tx", "exon", "cds"),
         type = c("tx", "exon", "cds"))
{
    gsubSQL <- function(long, short, sql) {
        gsub("LONG", long, gsub("SHORT", short, sql))
    }

    groupBy <- match.arg(groupBy)

    type <- match.arg(type)
    long <- ifelse(type == "tx", "transcript", type)
    short <- type

    ## create SQL query
    selectClause <-
      paste("SELECT DISTINCT SHORT_chrom, SHORT_start, SHORT_end,",
            "SHORT_strand, SHORT_name, LONG._SHORT_id AS SHORT_id")
    if (groupBy != "gene") {
        selectClause <-
          paste(selectClause, ", splicing._GROUPBY_id AS GROUPBY_id", sep = "")
    }
    if (type == "tx" || groupBy == "gene") {
        selectClause <- paste(selectClause, ", gene_id", sep = "")
    }
    fromClause <-
      paste("FROM LONG INNER JOIN LONG_rtree",
            "ON (LONG._SHORT_id=LONG_rtree._SHORT_id)")
    if (!(type == "tx" && groupBy == "gene")) {
        fromClause <-
          paste(fromClause,
                "INNER JOIN splicing",
                "ON (LONG._SHORT_id=splicing._SHORT_id)")
    }
    if (type == "tx" || groupBy == "gene") {
        selectClause <- paste(selectClause, ", gene_id", sep = "")
        fromClause <-
          paste(fromClause,
                " INNER JOIN gene",
                " ON (", ifelse(type == "tx", "transcript", "splicing"),
                "._tx_id=gene._tx_id)", sep = "")
    }
    whereClause <- "WHERE GROUPBY_id IS NOT NULL"
    orderByClause <-  "ORDER BY GROUPBY_id"
    if (type %in% c("exon", "cds") && groupBy == "tx") {
        orderByClause <- paste(orderByClause, ", splicing.exon_rank", sep = "")
    } else {
        orderByClause <-
          paste(orderByClause, ", LONG.SHORT_chrom, ",
                "LONG_rtree.SHORT_start, LONG_rtree.SHORT_end", sep = "")
    }

    sql <- paste(selectClause, fromClause, whereClause, orderByClause)
    sql <- gsub("LONG", long, sql)
    sql <- gsub("SHORT", short, sql)
    sql <- gsub("GROUPBY", groupBy, sql)

    ## get the data from the database
    ans <- dbGetQuery(txdb@conn, sql)

    ## create the GRanges object
    cols <- gsub("TYPE", type, c("TYPE_name", "TYPE_id"))
    if (type == "tx") {
        cols <- c(cols, "gene_id")
    }
    grngs <-
      GRanges(seqnames = ans[[paste(type, "_chrom", sep="")]],
              ranges = IRanges(start = ans[[paste(type, "_start", sep="")]],
                               end = ans[[paste(type, "_end", sep="")]]),
              strand = strand(ans[[paste(type, "_strand", sep="")]]),
              ans[cols])

    ## split by grouping variable
    split(grngs, ans[[paste(groupBy, "_id", sep="")]])
}


transcriptsByAlt <- function(txdb, groupBy = c("gene", "exon", "cds"))
{
    if(!is(txdb,"TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")
    groupBy <- match.arg(groupBy)
    .transcriptsORexonsORcdsBy(txdb, groupBy, "tx")
}

exonsByAlt <- function(txdb, groupBy = c("tx", "gene"))
{
    if(!is(txdb,"TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")
    groupBy <- match.arg(groupBy)
    .transcriptsORexonsORcdsBy(txdb, groupBy, "exon")
}

cdsByAlt <- function(txdb, groupBy = c("tx", "gene"))
{
    if(!is(txdb,"TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")
    groupBy <- match.arg(groupBy)
    .transcriptsORexonsORcdsBy(txdb, groupBy, "cds")
}
