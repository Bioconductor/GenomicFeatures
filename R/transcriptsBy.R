.transcriptsORexonsORcdsBy <-
function(txdb, by = c("gene", "tx", "exon", "cds"),
         type = c("tx", "exon", "cds"))
{
    if(!is(txdb,"TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")

    gsubSQL <- function(long, short, sql) {
        gsub("LONG", long, gsub("SHORT", short, sql))
    }

    by <- match.arg(by)

    type <- match.arg(type)
    long <- ifelse(type == "tx", "transcript", type)
    short <- type

    ## create SQL query
    selectClause <-
      paste("SELECT DISTINCT SHORT_chrom, SHORT_start, SHORT_end,",
            "SHORT_strand, SHORT_name, LONG._SHORT_id AS SHORT_id")
    if (by != "gene") {
        selectClause <-
          paste(selectClause, ", splicing._GROUPBY_id AS GROUPBY_id", sep = "")
    }
    if (type == "tx" || by == "gene") {
        selectClause <- paste(selectClause, ", gene_id", sep = "")
    }
    fromClause <-
      paste("FROM LONG INNER JOIN LONG_rtree",
            "ON (LONG._SHORT_id=LONG_rtree._SHORT_id)")
    if (!(type == "tx" && by == "gene")) {
        fromClause <-
          paste(fromClause,
                "INNER JOIN splicing",
                "ON (LONG._SHORT_id=splicing._SHORT_id)")
    }
    if (type == "tx" || by == "gene") {
        fromClause <-
          paste(fromClause,
                " INNER JOIN gene",
                " ON (", ifelse(type == "tx", "transcript", "splicing"),
                "._tx_id=gene._tx_id)", sep = "")
    }
    whereClause <- "WHERE GROUPBY_id IS NOT NULL"
    orderByClause <-  "ORDER BY GROUPBY_id"
    if (type %in% c("exon", "cds") && by == "tx") {
        orderByClause <- paste(orderByClause, ", exon_rank", sep = "")
    } else {
        orderByClause <-
          paste(orderByClause, ", SHORT_chrom, SHORT_strand, ",
                "SHORT_start, SHORT_end", sep = "")
    }

    sql <- paste(selectClause, fromClause, whereClause, orderByClause)
    sql <- gsub("LONG", long, sql)
    sql <- gsub("SHORT", short, sql)
    sql <- gsub("GROUPBY", by, sql)

    if (getOption("verbose", FALSE))
        cat("SQL QUERY: ", sql, "\n\n", sep = "")

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
    split(grngs, ans[[paste(by, "_id", sep="")]])
}


transcriptsBy <- function(txdb, by = c("gene", "exon", "cds"))
{
    .transcriptsORexonsORcdsBy(txdb, match.arg(by), "tx")
}

exonsBy <- function(txdb, by = c("tx", "gene"))
{
    .transcriptsORexonsORcdsBy(txdb, match.arg(by), "exon")
}

cdsBy <- function(txdb, by = c("tx", "gene"))
{
    .transcriptsORexonsORcdsBy(txdb, match.arg(by), "cds")
}
