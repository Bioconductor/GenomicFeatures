### =========================================================================
### Extraction of exonic and intronic parts
### -------------------------------------------------------------------------
###
### For all functions in this file, 'txdb' must be a TxDb object or any
### object that supports transcripts() and exonsBy() (e.g. EnsDb object).
###


### Return a GRanges object with 1 range per transcript and metadata columns
### tx_id, tx_name, and gene_id.
### If 'drop.geneless' is FALSE (the default) then the transcripts are returned
### in the same order as with transcripts(), which is expected to be by
### transcript id (tx_id). Otherwise they are ordered first by gene id
### (gene_id), then by transcript id.
.tidy_transcripts <- function(txdb, drop.geneless=FALSE)
{
    tx <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
    mcols(tx)$gene_id <- as.character(mcols(tx)$gene_id)
    if (drop.geneless) {
        gene_id <- mcols(tx)$gene_id
        tx_id <- mcols(tx)$tx_id
        tx <- tx[order(gene_id, tx_id, na.last=NA)]
    }
    tx
}

### Return a GRangesList object parallel to 'tx_ids'. The supplied 'tx_ids'
### must be a subset of 'mcols(transcripts(txdb))$tx_id'.
.exons_by_txids <- function(txdb, tx_ids)
{
    if (anyDuplicated(tx_ids))
        stop(wmsg("\"transcripts\" method for ", class(txdb), " objects ",
                  "seems broken, sorry"))
    ans <- exonsBy(txdb, by="tx")
    tx_ids <- as.character(tx_ids)
    ans_names <- names(ans)
    if (!identical(tx_ids, ans_names)) {
        m <- match(tx_ids, ans_names)
        if (anyNA(m))
            stop(wmsg("\"exonsBy\" method for ", class(txdb), " objects ",
                      "seems broken, sorry"))
        ans <- ans[m]
    }
    ans
}

### Return a GRanges object with 1 range per exon and metadata columns
### tx_id, tx_name, gene_id, exon_id, exon_name, and exon_rank.
### If 'drop.geneless' is FALSE (the default) then the exons are ordered first
### by transcript id (tx_id), then by exon rank (exon_rank). Otherwise they
### are ordered first by gene id (gene_id), then by transcript id, and then
### by exon rank.
.tidy_exons <- function(txdb, drop.geneless=FALSE)
{
    tx <- .tidy_transcripts(txdb, drop.geneless=drop.geneless)
    ex_by_tx <- .exons_by_txids(txdb, mcols(tx)$tx_id)

    ans <- unlist(ex_by_tx, use.names=FALSE)
    idx <- rep(seq_along(tx), lengths(ex_by_tx))
    mcols(ans) <- cbind(mcols(tx)[idx, , drop=FALSE], mcols(ans))
    ans
}

### Return a GRanges object with 1 range per intron and metadata columns
### tx_id, tx_name, and gene_id.
### If 'drop.geneless' is FALSE (the default) then the introns are ordered
### by transcript id (tx_id). Otherwise they are ordered first by gene id
### (gene_id), then by transcript id.
.tidy_introns <- function(txdb, drop.geneless=FALSE)
{
    tx <- .tidy_transcripts(txdb, drop.geneless=drop.geneless)
    ex_by_tx <- .exons_by_txids(txdb, mcols(tx)$tx_id)

    introns_by_tx <- psetdiff(tx, ex_by_tx)

    ans <- unlist(introns_by_tx, use.names=FALSE)
    idx <- rep(seq_along(tx), lengths(introns_by_tx))
    mcols(ans) <- mcols(tx)[idx, , drop=FALSE]
    ans
}

.break_in_parts <- function(x, linked.to.single.gene.only=FALSE)
{
    ans <- disjoin(x, with.revmap=TRUE)
    revmap <- mcols(ans)$revmap
    ans_mcols <- lapply(mcols(x),
                        function(col) unique(extractList(col, revmap)))
    mcols(ans) <- DataFrame(ans_mcols)
    if (linked.to.single.gene.only) {
        keep_idx <- which(elementNROWS(mcols(ans)$gene_id) == 1L)
        ans <- ans[keep_idx]
        mcols(ans)$gene_id <- as.character(mcols(ans)$gene_id)
    }
    ans
}

### Return a disjoint and strictly sorted GRanges object with 1 range per
### exonic part and metadata columns tx_id, tx_name, gene_id, exon_id,
### exon_name, and exon_rank.
exonicParts <- function(txdb, linked.to.single.gene.only=FALSE)
{
    if (!isTRUEorFALSE(linked.to.single.gene.only))
        stop("'linked.to.single.gene.only' must be TRUE or FALSE")
    ex <- .tidy_exons(txdb, drop.geneless=linked.to.single.gene.only)
    .break_in_parts(ex, linked.to.single.gene.only)
}

### Return a disjoint and strictly sorted GRanges object with 1 range per
### intronic part and metadata columns tx_id, tx_name, and gene_id.
intronicParts <- function(txdb, linked.to.single.gene.only=FALSE)
{
    if (!isTRUEorFALSE(linked.to.single.gene.only))
        stop("'linked.to.single.gene.only' must be TRUE or FALSE")
    introns <- .tidy_introns(txdb, drop.geneless=linked.to.single.gene.only)
    .break_in_parts(introns, linked.to.single.gene.only)
}

