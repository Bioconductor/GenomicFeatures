### =========================================================================
### transcriptLengths()
### -------------------------------------------------------------------------


.match_and_check <- function(rglist_names, tx_id)
{
    if (is.null(rglist_names))
        stop(wmsg("internal error in transcriptLengths(): ",
                  "no names on 'rglist'"))
    m <- match(rglist_names, tx_id)
    if (any(is.na(m)))
        stop(wmsg("internal error in transcriptLengths(): ",
                  "some 'rglist' names cannot be mapped to 'tx_id'"))
    m
}

### 'rglist' must be a named RangesList or GRangesList.
### 'tx_id' must be a character vector.
.eltNROWS <- function(rglist, tx_id)
{
    ans <- integer(length(tx_id))
    m <- .match_and_check(names(rglist), tx_id)
    ans[m] <- elementNROWS(rglist)
    ans
}

.sum_width <- function(rglist, tx_id)
{
    ans <- integer(length(tx_id))
    m <- .match_and_check(names(rglist), tx_id)
    ans[m] <- sum(width(rglist))
    ans
}

### The returned data frame has 1 row per transcript returned by
### 'transcripts(txdb)' and in the same order.
### NOTES:
### - The functions only accepts a TxDb object for now. We'll make it
###   a generic function when we need to support other types of input.
### - The function could probably be made much faster by querying the
###   TxDb object directly in SQL instead of calling exonsBy(), cdsBy(),
###   fiveUTRsByTranscript(), and threeUTRsByTranscript() successively.
transcriptLengths <- function(txdb, with.cds_len=FALSE,
                                    with.utr5_len=FALSE, with.utr3_len=FALSE)
{
    if (!isTRUEorFALSE(with.cds_len))
        stop("'with.cds_len' must be TRUE or FALSE")
    if (!isTRUEorFALSE(with.utr5_len))
        stop("'with.utr5_len' must be TRUE or FALSE")
    if (!isTRUEorFALSE(with.cds_len))
        stop("'with.utr3_len' must be TRUE or FALSE")
    tx <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
    ans <- mcols(tx)
    ans$gene_id <- as.character(ans$gene_id)
    tx_id <- as.character(ans$tx_id)  # because match() will want a character

    rg_by_tx <- exonsBy(txdb, by="tx")
    ans$nexon <- .eltNROWS(rg_by_tx, tx_id)
    ans$tx_len <- .sum_width(rg_by_tx, tx_id)
    if (with.cds_len) {
        rg_by_tx <- cdsBy(txdb, by="tx")
        ans$cds_len <- .sum_width(rg_by_tx, tx_id)
    }
    if (with.utr5_len) {
        rg_by_tx <- fiveUTRsByTranscript(txdb)
        ans$utr5_len <- .sum_width(rg_by_tx, tx_id)
    }
    if (with.utr3_len) {
        rg_by_tx <- threeUTRsByTranscript(txdb)
        ans$utr3_len <- .sum_width(rg_by_tx, tx_id)
    }
    as.data.frame(ans)
}

