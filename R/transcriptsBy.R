### =========================================================================
### The "features grouped by" extractors
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### id2name().
###

id2name <- function(txdb, feature.type=c("tx", "exon", "cds"))
{
    if(!is(txdb,"TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")
    feature.type <- match.arg(feature.type)
    tablename <- switch(feature.type, tx="transcript", exon="exon", cds="cds")
    what_cols <- paste(c("_", ""), feature.type, c("_id", "_name"), sep="")
    SQL <- paste("SELECT",
                 paste(what_cols, collapse=", "),
                 "FROM", tablename)
    data <- dbEasyQuery(txdbConn(txdb), SQL)
    ans <- data[[what_cols[2L]]]
    names(ans) <- as.character(data[[what_cols[1L]]])
    ans
}

.set.group.names <- function(ans, use.names, txdb, by)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (!use.names)
        return(ans)
    if (by == "gene")
        stop("'use.names=TRUE' cannot be used when grouping by gene")
    names(ans) <- id2name(txdb, feature.type=by)[names(ans)]
    if (any(is.na(names(ans))) || any(duplicated(names(ans))))
        warning("some group names are NAs or duplicated")
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### transcriptsBy(), exonsBy() and cdsBy().
###

.featuresBy <- function(txdb, by, type,
                        distinct=FALSE,
                        splicing_in_join=TRUE,
                        gene_in_join=FALSE,
                        order_by_exon_rank=TRUE,
                        use.names=FALSE)
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

    ## get the data from the database
    data <- dbEasyQuery(txdbConn(txdb), sql)

    ## create the GRanges object
    cols <- gsub("TYPE", type, c("TYPE_id", "TYPE_name"))
    if (order_by_exon_rank)
        cols <- c(cols, "exon_rank")
    seqinfo <- seqinfo(txdb)
    grngs <-
      GRanges(seqnames =
              factor(data[[paste(type, "_chrom", sep="")]],
                     levels = seqlevels(seqinfo)),
              ranges = IRanges(start = data[[paste(type, "_start", sep="")]],
                               end = data[[paste(type, "_end", sep="")]]),
              strand = strand(data[[paste(type, "_strand", sep="")]]),
              data[cols])
    seqinfo(grngs) <- seqinfo

    ## split by grouping variable
    ans <- split(grngs, data[[paste(by, "_id", sep="")]])
    .set.group.names(ans, use.names, txdb, by)
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

setMethod("transcriptsBy", "TranscriptDb",
    function(x, by=c("gene", "exon", "cds"), use.names=FALSE)
    {
        by <- match.arg(by)
        distinct <- splicing_in_join <- by != "gene"
        gene_in_join <- by == "gene"
        .featuresBy(x, by, "tx",
                    distinct=distinct,
                    splicing_in_join=splicing_in_join,
                    gene_in_join=gene_in_join,
                    order_by_exon_rank=FALSE,
                    use.names=use.names)
    }
)

setMethod("exonsBy", "TranscriptDb",
    function(x, by=c("tx", "gene"), use.names=FALSE)
    {
        by <- match.arg(by)
        distinct <- gene_in_join <- by == "gene"
        order_by_exon_rank <- by == "tx"
        .featuresBy(x, by, "exon",
                    distinct=distinct,
                    splicing_in_join=TRUE,
                    gene_in_join=gene_in_join,
                    order_by_exon_rank=order_by_exon_rank,
                    use.names=use.names)
    }
)

setMethod("cdsBy", "TranscriptDb",
    function(x, by=c("tx", "gene"), use.names=FALSE)
    {
        by <- match.arg(by)
        distinct <- gene_in_join <- by == "gene"
        order_by_exon_rank <- by == "tx"
        .featuresBy(x, by, "cds",
                    distinct=distinct,
                    splicing_in_join=TRUE,
                    gene_in_join=gene_in_join,
                    order_by_exon_rank=order_by_exon_rank,
                    use.names=use.names)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### intronsByTranscript().
###

setMethod("intronsByTranscript", "TranscriptDb",
    function(x, use.names=FALSE)
    {
        tx <- transcripts(x)
        exn <- exonsBy(x)
        tx <- tx[match(names(exn), elementMetadata(tx)[,"tx_id"])]
        ans <- psetdiff(tx, exn)
        .set.group.names(ans, use.names, x, "tx")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fiveUTRsByTranscript() and threeUTRsByTranscript().
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
    ans <- dbEasyQuery(txdbConn(txdb), sql)
    if (translated.transcripts.only) {
        ids <- unique(ans$tx_id[!is.na(ans$cds_id)])
        ans <- ans[ans$tx_id %in% ids, ]
    }
    ans
}

.makeUTRsByTranscript <- function(splicings, utr_start, utr_end, seqinfo)
{
    seqlevels <- seqlevels(seqinfo)
    grg <- GRanges(seqnames=factor(splicings$exon_chrom, levels=seqlevels),
                   ranges=IRanges(start=utr_start, end=utr_end),
                   strand=strand(splicings$exon_strand),
                   exon_id=splicings$exon_id,
                   exon_name=splicings$exon_name,
                   exon_rank=splicings$exon_rank)
    seqinfo(grg) <- seqinfo
    idx <- width(grg) != 0L  # drop 0-width UTRs
    split(grg[idx], splicings$tx_id[idx])
}

setMethod("fiveUTRsByTranscript", "TranscriptDb",
    function(x, use.names=FALSE)
    {
        splicings <- .getFullSplicings(x, translated.transcripts.only=TRUE)

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
        seqinfo <- seqinfo(x)
        ans <- .makeUTRsByTranscript(splicings, utr_start, utr_end, seqinfo)
        .set.group.names(ans, use.names, x, "tx")
    }
)

setMethod("threeUTRsByTranscript", "TranscriptDb",
    function(x, use.names=FALSE)
    {
        splicings <- .getFullSplicings(x, translated.transcripts.only=TRUE)

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
        seqinfo <- seqinfo(x)
        ans <- .makeUTRsByTranscript(splicings, utr_start, utr_end, seqinfo)
        .set.group.names(ans, use.names, x, "tx")
    }
)
