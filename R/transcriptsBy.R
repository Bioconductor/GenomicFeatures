### =========================================================================
### The "features grouped by" extractors
### -------------------------------------------------------------------------


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

## helper to translate back to what is expected from seqinfo()
.translateToSeqInfo <- function(txdb, x){
    tr <- load_chrominfo(txdb, set.col.class=TRUE)$chrom[txdb$user2seqlevels0]
    names(tr) <- txdb$user_seqlevels
    idx <- match(x, tr)
    names(tr)[idx]
}

## !!!
## TODO: I think that this helper is screwing things up
.baseNamedActiveSeqs <- function(txdb){
    trueNames <- load_chrominfo(txdb, set.col.class=TRUE)$chrom
    actSqs <- .isActiveSeq(txdb)
    names(actSqs) <- trueNames[txdb$user2seqlevels0] ## limit result to these.
    actSqs
}

.getOnlyActiveSeqs <- function(txdb){
    actSqs <- .baseNamedActiveSeqs(txdb)
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
    if(!is(txdb, "TxDb"))
        stop("'txdb' must be a TxDb object")

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
    whereSeqsClause <- paste0("AND SHORT_chrom IN ('",
                              paste(.getOnlyActiveSeqs(txdb), collapse="','")
                              ,"')")
    
    sql <- paste(selectClause, fromClause, whereClause, whereSeqsClause,
                 orderByClause)
    sql <- gsub("LONG", long, sql)
    sql <- gsub("SHORT", short, sql)
    sql <- gsub("GROUPBY", by, sql)

    ## get the data from the database
    data <- queryAnnotationDb(txdb, sql)
    ## seqnames may be out of sync with expected results. Massage back.
    chromName <- paste0(type,"_chrom")
    data[[chromName]] <- .translateToSeqInfo(txdb, data[[chromName]])

    ## create the GRanges object
    cols <- paste0(type, c("_id", "_name"))
    if (order_by_exon_rank)
        cols <- c(cols, "exon_rank")
    activeNames <- names(.isActiveSeq(txdb))[.isActiveSeq(txdb)]
    seqinfo <- seqinfo(txdb)[activeNames]
    grngs <- GRanges(seqnames = factor(
                       data[[paste0(type, "_chrom")]],
                       levels = activeNames),
                     ranges = IRanges(
                       start = data[[paste0(type, "_start")]],
                       end = data[[paste0(type, "_end")]]),
                     strand = strand(data[[paste0(type, "_strand")]]),
                     data[cols],
                     seqinfo = seqinfo)

    ## split by grouping variable
    ans <- split(grngs, data[[paste0(by, "_id")]])
    ans <- .set.group.names(ans, use.names, txdb, by)
    .assignMetadataList(ans, txdb)
}

setGeneric("transcriptsBy", signature="x",
    function(x, by=c("gene", "exon", "cds"), ...)
        standardGeneric("transcriptsBy")
)

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

setMethod("transcriptsBy", "TxDb",
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

setGeneric("exonsBy", signature="x",
    function(x, by=c("tx", "gene"), ...) standardGeneric("exonsBy")
)

setMethod("exonsBy", "TxDb",
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

setGeneric("cdsBy", signature="x",
    function(x, by=c("tx", "gene"), ...) standardGeneric("cdsBy")
)

setMethod("cdsBy", "TxDb",
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

setGeneric("intronsByTranscript",
    function(x, ...) standardGeneric("intronsByTranscript")
)

setMethod("intronsByTranscript", "TxDb",
    function(x, use.names=FALSE)
    {
        tx <- transcripts(x)
        exn <- exonsBy(x)
        tx <- tx[match(names(exn), mcols(tx)[,"tx_id"])]
        ans <- psetdiff(tx, exn)
        .set.group.names(ans, use.names, x, "tx")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fiveUTRsByTranscript() and threeUTRsByTranscript().
###

.getSplicingsForTranscriptsWithCDSs <- function(txdb)
{
    ans <- load_splicings(txdb)
    ids <- unique(ans$tx_id[!is.na(ans$cds_id)])
    ans[ans$tx_id %in% ids, ]
}

### 'tx_id': character or integer vector with runs of identical elements (one
### run per transcript, and, within each run, one element per exon).
### 'exons_with_cds': integer vector containing the indices of the elements
### in 'tx_id' corresponding to exons with a CDS. The indices must be sorted
### in ascending order.
### Returns for each transcript the indices of the first exon with a CDS, plus
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
    S4Vectors:::fancy_mseq(lengths, offset=offset)
}

### 'tx_id', 'exons_with_cds': same as for .exons_with_5utr().
### Returns for each transcript the indices of the last exon with a CDS, plus
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
    S4Vectors:::fancy_mseq(lengths, offset=offset)
}

.makeUTRsByTranscript <- function(txdb, splicings, utr_start, utr_end)
{
    seqinfo0 <- get_TxDb_seqinfo0(txdb)
    cols <- paste0("exon_", c("id", "name", "rank"))
    gr <- GRanges(seqnames = factor(
                     splicings$exon_chrom,
                     levels =seqlevels(seqinfo0)),
                  ranges = IRanges(start = utr_start, end = utr_end),
                  strand = strand(splicings$exon_strand),
                  splicings[cols],
                  seqinfo = seqinfo0)
    idx <- width(gr) != 0L  # drop 0-width UTRs
    grl <- split(gr[idx], splicings$tx_id[idx])
    keep_user_seqlevels_from_TxDb(grl, txdb)
}

.make5UTRsByTranscript <- function(txdb, splicings, use.names=FALSE)
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
    ans <- .makeUTRsByTranscript(txdb, splicings, utr_start, utr_end)
    ans <- .set.group.names(ans, use.names, txdb, "tx")
    ans <- .assignMetadataList(ans, txdb)
    ans
}

.make3UTRsByTranscript <- function(txdb, splicings, use.names=FALSE)
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
    ans <- .makeUTRsByTranscript(txdb, splicings, utr_start, utr_end)
    ans <- .set.group.names(ans, use.names, txdb, "tx")            
    ans <- .assignMetadataList(ans, txdb)
    ans
}

setGeneric("fiveUTRsByTranscript",
    function(x, ...) standardGeneric("fiveUTRsByTranscript")
)

setMethod("fiveUTRsByTranscript", "TxDb",
    function(x, use.names=FALSE)
    {
        splicings <- .getSplicingsForTranscriptsWithCDSs(x)
        .make5UTRsByTranscript(x, splicings, use.names=use.names)
    }
)

setGeneric("threeUTRsByTranscript",
    function(x, ...) standardGeneric("threeUTRsByTranscript")
)

setMethod("threeUTRsByTranscript", "TxDb",
    function(x, use.names=FALSE)
    {
        splicings <- .getSplicingsForTranscriptsWithCDSs(x)
        .make3UTRsByTranscript(x, splicings, use.names=use.names)
    }
)

