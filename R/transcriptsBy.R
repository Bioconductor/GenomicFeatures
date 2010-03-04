.featuresBy <- function(txdb, by, type,
                        distinct=FALSE,
                        splicing_in_join=TRUE,
                        gene_in_join=FALSE,
                        order_by_exon_rank=TRUE)
{
    if(!is(txdb,"TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")

    long <- ifelse(type == "tx", "transcript", type)
    short <- type

    ## create SQL query
    selectClause <- "SELECT"
    if (distinct)
        selectClause <- paste(selectClause, "DISTINCT")
    selectClause <-
      paste(selectClause, "SHORT_chrom, SHORT_start, SHORT_end,",
            "SHORT_strand, SHORT_name, LONG._SHORT_id AS SHORT_id")
    if (by == "gene") {
        selectClause <- paste(selectClause, ", gene_id", sep = "")
    } else {
        selectClause <-
          paste(selectClause, ", splicing._GROUPBY_id AS GROUPBY_id", sep = "")
    }
    fromClause <- "FROM LONG"
    if (splicing_in_join) {
        fromClause <-
          paste(fromClause,
                "INNER JOIN splicing",
                "ON (LONG._SHORT_id=splicing._SHORT_id)")
    }
    if (gene_in_join) {
        fromClause <-
          paste(fromClause,
                " INNER JOIN gene",
                " ON (", ifelse(type == "tx", "transcript", "splicing"),
                "._tx_id=gene._tx_id)", sep = "")
    }
    whereClause <- "WHERE GROUPBY_id IS NOT NULL"
    orderByClause <-  "ORDER BY GROUPBY_id"
    if (order_by_exon_rank) {
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
    grngs <-
      GRanges(seqnames = ans[[paste(type, "_chrom", sep="")]],
              ranges = IRanges(start = ans[[paste(type, "_start", sep="")]],
                               end = ans[[paste(type, "_end", sep="")]]),
              strand = strand(ans[[paste(type, "_strand", sep="")]]),
              ans[cols])

    ## split by grouping variable
    split(grngs, ans[[paste(by, "_id", sep="")]])
}

###                    use  splicing      gene
###   type    by  DISTINCT   in JOIN   in JOIN   ORDER BY
##    ----  ----  --------  --------  --------  ---------
###     tx  gene        no        no       yes      locus
###     tx  exon        no       yes        no      locus
###     tx   cds        no       yes        no      locus
###   exon    tx        no       yes        no  exon_rank
###   exon  gene       yes       yes       yes      locus
###    cds    tx        no       yes        no  exon_rank
###    cds  gene       yes       yes       yes      locus

transcriptsBy <- function(txdb, by = c("gene", "exon", "cds"))
{
    by <- match.arg(by)
    splicing_in_join <- by != "gene"
    gene_in_join <- by == "gene"
    .featuresBy(txdb, by, "tx",
                distinct=FALSE,
                splicing_in_join=splicing_in_join,
                gene_in_join=gene_in_join,
                order_by_exon_rank=FALSE)
}

exonsBy <- function(txdb, by = c("tx", "gene"))
{
    by <- match.arg(by)
    distinct <- gene_in_join <- by == "gene"
    order_by_exon_rank <- by == "tx"
    .featuresBy(txdb, by, "exon",
                distinct=distinct,
                splicing_in_join=TRUE,
                gene_in_join=gene_in_join,
                order_by_exon_rank=order_by_exon_rank)
}

cdsBy <- function(txdb, by = c("tx", "gene"))
{
    by <- match.arg(by)
    distinct <- gene_in_join <- by == "gene"
    order_by_exon_rank <- by == "tx"
    .featuresBy(txdb, by, "cds",
                distinct=distinct,
                splicing_in_join=TRUE,
                gene_in_join=gene_in_join,
                order_by_exon_rank=order_by_exon_rank)
}

