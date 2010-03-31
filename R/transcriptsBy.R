### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### transcriptsBy(), exonsBy() and cdsBy().
###

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
      paste(selectClause, "LONG._SHORT_id AS SHORT_id, SHORT_name,",
            "SHORT_chrom, SHORT_strand, SHORT_start, SHORT_end")
    if (by == "gene") {
        selectClause <- paste(selectClause, ", gene_id", sep = "")
    } else if (splicing_in_join) {
        selectClause <-
          paste(selectClause, ", splicing._GROUPBY_id AS GROUPBY_id", sep = "")
        if (order_by_exon_rank)
            selectClause <- paste(selectClause, ", exon_rank", sep = "")
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
    ans <- dbGetQuery(txdbConn(txdb), sql)

    ## create the GRanges object
    cols <- gsub("TYPE", type, c("TYPE_id", "TYPE_name"))
    if (order_by_exon_rank)
        cols <- c(cols, "exon_rank")
    seqlengths <- seqlengths(txdb)
    grngs <-
      GRanges(seqnames =
              factor(ans[[paste(type, "_chrom", sep="")]],
                     levels = names(seqlengths)),
              ranges = IRanges(start = ans[[paste(type, "_start", sep="")]],
                               end = ans[[paste(type, "_end", sep="")]]),
              strand = strand(ans[[paste(type, "_strand", sep="")]]),
              ans[cols],
              seqlengths = seqlengths)

    ## split by grouping variable
    split(grngs, ans[[paste(by, "_id", sep="")]])
}

###                    use  splicing      gene
###   type    by  DISTINCT   in JOIN   in JOIN   ORDER BY
##    ----  ----  --------  --------  --------  ---------
###     tx  gene        no        no       yes      locus
###     tx  exon       yes       yes        no      locus
###     tx   cds       yes       yes        no      locus
###   exon    tx        no       yes        no  exon_rank
###   exon  gene       yes       yes       yes      locus
###    cds    tx        no       yes        no  exon_rank
###    cds  gene       yes       yes       yes      locus

transcriptsBy <- function(txdb, by = c("gene", "exon", "cds"))
{
    by <- match.arg(by)
    distinct <- splicing_in_join <- by != "gene"
    gene_in_join <- by == "gene"
    .featuresBy(txdb, by, "tx",
                distinct=distinct,
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fiveUTRsByTranscripts() and threeUTRsByTranscripts().
###

.getFullSplicings <- function(txdb, translated.transcripts.only=FALSE)
{
    ORDER_BY <- "ORDER BY tx_chrom, tx_strand, tx_start, tx_end, tx_id"
    sql <- paste(
        "SELECT transcript._tx_id AS tx_id, exon_rank,",
        " splicing._exon_id AS exon_id,",
        " exon_name, exon_chrom, exon_strand, exon_start, exon_end,",
        " splicing._cds_id AS cds_id,",
        " cds_start, cds_end",
        "FROM transcript",
        " INNER JOIN splicing",
        "  ON (tx_id=splicing._tx_id)",
        " INNER JOIN exon",
        "  ON (exon_id=exon._exon_id)",
        " LEFT JOIN cds",
        "  ON (cds_id=cds._cds_id)",
        ORDER_BY, ", exon_rank")
    ans <- dbGetQuery(txdbConn(txdb), sql)
    if (translated.transcripts.only) {
        ids <- unique(ans$tx_id[!is.na(ans$cds_id)])
        ans <- ans[ans$tx_id %in% ids, ]
    }
    ans
}

.makeUTRsByTranscripts <- function(splicings, utr_start, utr_end, seqlengths)
{
    seqlevels <- names(seqlengths)
    grg <- GRanges(seqnames=factor(splicings$exon_chrom, levels=seqlevels),
                   ranges=IRanges(start=utr_start, end=utr_end),
                   strand=strand(splicings$exon_strand),
                   exon_id=splicings$exon_id,
                   exon_name=splicings$exon_name,
                   exon_rank=splicings$exon_rank,
                   seqlengths=seqlengths)
    idx <- width(grg) != 0L  # drop 0-width UTRs
    split(grg[idx], splicings$tx_id[idx])
}

fiveUTRsByTranscripts <- function(txdb)
{
    splicings <- .getFullSplicings(txdb, translated.transcripts.only=TRUE)

    ## For each transcript, we keep only the first row with a CDS plus all
    ## previous rows (if any).
    if (nrow(splicings) != 0L) {
        cdslist <- split(splicings$cds_id, splicings$tx_id)
        tmp <- lapply(cdslist,
                 function(cds_id)
                 {
                     W <- which(!is.na(cds_id))
                     L <- W[1L]
                     rep.int(c(TRUE, FALSE), c(L, length(cds_id)-L))
                 })
        idx <- unsplit(tmp, splicings$tx_id)
        splicings <- splicings[idx, ]
    }

    ## Compute the UTR starts/ends.
    utr_start <- splicings$exon_start
    utr_end <- splicings$exon_end
    idx1 <- !is.na(splicings$cds_id)
    idx <- idx1 & (splicings$exon_strand == "+")
    utr_end[idx] <- splicings$cds_start[idx] - 1L
    idx <- idx1 & (splicings$exon_strand == "-")
    utr_start[idx] <- splicings$cds_end[idx] + 1L

    ## Make and return the GRangesList object.
    seqlengths <- seqlengths(txdb)
    .makeUTRsByTranscripts(splicings, utr_start, utr_end, seqlengths)
}

threeUTRsByTranscripts <- function(txdb)
{
    splicings <- .getFullSplicings(txdb, translated.transcripts.only=TRUE)

    ## For each transcript, we keep only the last row with a CDS plus all
    ## following rows (if any).
    if (nrow(splicings) != 0L) {
        cdslist <- split(splicings$cds_id, splicings$tx_id)
        tmp <- lapply(cdslist,
                 function(cds_id)
                 {
                     W <- which(!is.na(cds_id))
                     L <- W[length(W)]
                     rep.int(c(FALSE, TRUE), c(L-1L, length(cds_id)-L+1L))
                 })
        idx <- unsplit(tmp, splicings$tx_id)
        splicings <- splicings[idx, ]
    }

    ## Compute the UTR starts/ends.
    utr_start <- splicings$exon_start
    utr_end <- splicings$exon_end
    idx1 <- !is.na(splicings$cds_id)
    idx <- idx1 & (splicings$exon_strand == "+")
    utr_start[idx] <- splicings$cds_end[idx] + 1L
    idx <- idx1 & (splicings$exon_strand == "-")
    utr_end[idx] <- splicings$cds_start[idx] - 1L

    ## Make and return the GRangesList object.
    seqlengths <- seqlengths(txdb)
    .makeUTRsByTranscripts(splicings, utr_start, utr_end, seqlengths)
}

