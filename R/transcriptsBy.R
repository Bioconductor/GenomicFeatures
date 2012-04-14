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
    data <- dbEasyQuery(AnnotationDbi:::dbConn(txdb), SQL)
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

.getOnlyActiveSeqs <- function(txdb){
    actSqs <- isActiveSeq(txdb)
    names(actSqs)[actSqs]
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
    orderByClause <- "ORDER BY GROUPBY_id"
    if (order_by_exon_rank) {
        orderByClause <- paste(orderByClause, ", exon_rank", sep = "")
    } else {
        orderByClause <-
          paste(orderByClause, ", SHORT_chrom, SHORT_strand, ",
                "SHORT_start, SHORT_end", sep = "")
    }
    whereSeqsClause <- paste("AND SHORT_chrom IN ('",
                             paste(.getOnlyActiveSeqs(txdb),collapse="','")
                             ,"')", sep="")
    
    sql <- paste(selectClause, fromClause, whereClause, whereSeqsClause,
                 orderByClause)
    sql <- gsub("LONG", long, sql)
    sql <- gsub("SHORT", short, sql)
    sql <- gsub("GROUPBY", by, sql)

    ## get the data from the database
    data <- dbEasyQuery(AnnotationDbi:::dbConn(txdb), sql)

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
    ## Filter seqinfo
    isActSeq <- isActiveSeq(txdb)
    seqinfo(grngs) <- seqinfo
    seqlevels(grngs) <- names(isActSeq)[isActSeq]

    ## split by grouping variable
    ans <- split(grngs, data[[paste(by, "_id", sep="")]])
    ans <- .set.group.names(ans, use.names, txdb, by)
    ans <- .assignMetadataList(ans, txdb)
    ans
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

.getSplicingsForTranscriptsWithCDSs <- function(txdb)
{
    ans <- getSplicings(txdb)
    ids <- unique(ans$tx_id[!is.na(ans$cds_id)])
    ans <- ans[ans$tx_id %in% ids, ]

    ## modify results to not respect our activeSeqs mask.
    isActSeq <- isActiveSeq(txdb)
    ## remove unwanted stuff from df 
    remove <- names(isActSeq)[isActSeq==FALSE]
    ans[!(ans$exon_chrom %in% remove),]
}

### 'tx_id': character or integer vector with runs of identical elements (one
### run per transcript, and, within each run, one element per exon).
### 'exons_with_cds': integer vector containing the indices of the elements
### in 'tx_id' corresponding to exons with a CDS. The indices must be sorted
### in ascending order.
### Returns the indices of the first exon with a CDS for every transcript, plus
### the indices of the preceding exons in the transcript.
### Note that, for each transcript, exons that have a 5' UTR are all the exons
### before the first exon with a CDS, including the first exon with a CDS (even
### though this one might actually have a 0-width 5' UTR but we will take care
### of this later).
.exons_with_5utr <- function(tx_id, exons_with_cds)
{
    tx_id_rle <- Rle(tx_id)
    nexon <- runLength(tx_id_rle)
    ntx <- length(nexon)
    pseudo_tx_id <- rep.int(seq_len(ntx), nexon)
    first_exon_with_cds <- integer(ntx)
    exons_with_cds <- rev(exons_with_cds)
    first_exon_with_cds[pseudo_tx_id[exons_with_cds]] <- exons_with_cds
    offset <- cumsum(c(0L, nexon[-length(nexon)]))
    lengths <- first_exon_with_cds - offset
    IRanges:::fancy_mseq(lengths, offset=offset)
}

### 'tx_id', 'exons_with_cds': same as for .exons_with_5utr().
### Returns the indices of the last exon with a CDS for every transcript, plus
### the indices of the following exons in the transcript.
### Note that, for each transcript, exons that have a 3' UTR are all the exons
### after the last exon with a CDS, including the last exon with a CDS (even
### though this one might actually have a 0-width 3' UTR but we will take care
### of this later).
.exons_with_3utr <- function(tx_id, exons_with_cds)
{
    tx_id_rle <- Rle(tx_id)
    nexon <- runLength(tx_id_rle)
    ntx <- length(nexon)
    pseudo_tx_id <- rep.int(seq_len(ntx), nexon)
    last_exon_with_cds <- integer(ntx)
    last_exon_with_cds[pseudo_tx_id[exons_with_cds]] <- exons_with_cds
    offset <- last_exon_with_cds - 1L
    lengths <- cumsum(nexon) - offset
    IRanges:::fancy_mseq(lengths, offset=offset)
}

.makeUTRsByTranscript <- function(x, splicings, utr_start, utr_end)
{
    seqinfo <- seqinfo(x)
    grg <- GRanges(seqnames=factor(splicings$exon_chrom, 
                                   levels=seqlevels(seqinfo)),
                   ranges=IRanges(start=utr_start, end=utr_end),
                   strand=strand(splicings$exon_strand),
                   exon_id=splicings$exon_id,
                   exon_name=splicings$exon_name,
                   exon_rank=splicings$exon_rank)
    ## Then clean up the seqinfo
    isActSeq <- isActiveSeq(x)
    seqinfo(grg) <- seqinfo
    seqlevels(grg) <- names(isActSeq)[isActSeq]
    idx <- width(grg) != 0L  # drop 0-width UTRs
    split(grg[idx], splicings$tx_id[idx])
}

.make5UTRsByTranscript <- function(x, splicings, use.names=FALSE)
{
    exons_with_cds <- which(!is.na(splicings$cds_id))
    idx <- .exons_with_5utr(splicings$tx_id, exons_with_cds)
    splicings <- splicings[idx, ]

    ## Compute the UTR starts/ends.
    utr_start <- splicings$exon_start
    utr_end <- splicings$exon_end
    idx1 <- !is.na(splicings$cds_id)
    idx <- idx1 & (splicings$exon_strand == "+")
    utr_end[idx] <- splicings$cds_start[idx] - 1L
    idx <- idx1 & (splicings$exon_strand == "-")
    utr_start[idx] <- splicings$cds_end[idx] + 1L

    ## split by grouping variable
    ans <- .makeUTRsByTranscript(x, splicings, utr_start, utr_end)
    ans <- .set.group.names(ans, use.names, x, "tx")
    ans <- .assignMetadataList(ans, x)
    ans
}

.make3UTRsByTranscript <- function(x, splicings, use.names=FALSE)
{
    exons_with_cds <- which(!is.na(splicings$cds_id))
    idx <- .exons_with_3utr(splicings$tx_id, exons_with_cds)
    splicings <- splicings[idx, ]

    ## Compute the UTR starts/ends.
    utr_start <- splicings$exon_start
    utr_end <- splicings$exon_end
    idx1 <- !is.na(splicings$cds_id)
    idx <- idx1 & (splicings$exon_strand == "+")
    utr_start[idx] <- splicings$cds_end[idx] + 1L
    idx <- idx1 & (splicings$exon_strand == "-")
    utr_end[idx] <- splicings$cds_start[idx] - 1L

    ## split by grouping variable
    ans <- .makeUTRsByTranscript(x, splicings, utr_start, utr_end)
    ans <- .set.group.names(ans, use.names, x, "tx")            
    ans <- .assignMetadataList(ans, x)
    ans
}

setMethod("fiveUTRsByTranscript", "TranscriptDb",
    function(x, use.names=FALSE)
    {
        splicings <- .getSplicingsForTranscriptsWithCDSs(x)
        .make5UTRsByTranscript(x, splicings, use.names=use.names)
    }
)

setMethod("threeUTRsByTranscript", "TranscriptDb",
    function(x, use.names=FALSE)
    {
        splicings <- .getSplicingsForTranscriptsWithCDSs(x)
        .make3UTRsByTranscript(x, splicings, use.names=use.names)
    }
)

