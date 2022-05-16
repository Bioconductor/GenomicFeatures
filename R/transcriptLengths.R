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

### 'rglist' must be a named IntegerRangesList or GRangesList.
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


#' Extract the transcript lengths (and other metrics) from a TxDb object
#' 
#' The \code{transcriptLengths} function extracts the transcript lengths from a
#' \link{TxDb} object. It also returns the CDS and UTR lengths for each
#' transcript if the user requests them.
#' 
#' All the lengths are counted in number of nucleotides.
#' 
#' The length of a processed transcript is just the sum of the lengths of its
#' exons. This should not be confounded with the length of the stretch of DNA
#' transcribed into RNA (a.k.a. transcription unit), which can be obtained with
#' \code{width(transcripts(txdb))}.
#' 
#' @param txdb A \link{TxDb} object.
#' @param with.cds_len,with.utr5_len,with.utr3_len \code{TRUE} or \code{FALSE}.
#' Whether or not to also extract and return the CDS, 5' UTR, and 3' UTR
#' lengths for each transcript.
#' @param \dots Additional arguments used by \code{transcripts} and other
#' accessor functions.
#' @return A data frame with 1 row per transcript. The rows are guaranteed to
#' be in the same order as the elements of the \link[GenomicRanges]{GRanges}
#' object returned by \code{\link{transcripts}(txdb)}.  The data frame has
#' between 5 and 8 columns, depending on what the user requested via the
#' \code{with.cds_len}, \code{with.utr5_len}, and \code{with.utr3_len}
#' arguments.
#' 
#' The first 3 columns are the same as the metadata columns of the object
#' returned by \preformatted{ transcripts(txdb, columns=c("tx_id", "tx_name",
#' "gene_id")) } that is: \itemize{ \item \code{tx_id}: The internal transcript
#' ID. This ID is unique within the scope of the \link{TxDb} object. It is not
#' an official or public ID (like an Ensembl or FlyBase ID) or an Accession
#' number, so it cannot be used to lookup the transcript in public data bases
#' or in other \link{TxDb} objects. Furthermore, this ID could change when
#' re-running the code that was used to make the \link{TxDb} object.  \item
#' \code{tx_name}: An official/public transcript name or ID that can be used to
#' lookup the transcript in public data bases or in other \link{TxDb} objects.
#' This column is not guaranteed to contain unique values and it can contain
#' NAs.  \item \code{gene_id}: The official/public ID of the gene that the
#' transcript belongs to. Can be NA if the gene is unknown or if the transcript
#' is not considered to belong to a gene.  }
#' 
#' The other columns are quantitative: \itemize{ \item \code{nexon}: The number
#' of exons in the transcript.  \item \code{tx_len}: The length of the
#' processed transcript.  \item \code{cds_len}: [optional] The length of the
#' CDS region of the processed transcript.  \item \code{utr5_len}: [optional]
#' The length of the 5' UTR region of the processed transcript.  \item
#' \code{utr3_len}: [optional] The length of the 3' UTR region of the processed
#' transcript.  }
#' @author Hervé Pagès
#' @seealso \itemize{ \item \code{\link{transcripts}},
#' \code{\link{transcriptsBy}}, and \code{\link{transcriptsByOverlaps}}, for
#' extracting genomic feature locations from a \link{TxDb}-like object.
#' 
#' \item \code{\link{exonicParts}} and \code{\link{intronicParts}} for
#' extracting non-overlapping exonic or intronic parts from a TxDb-like object.
#' 
#' \item \code{\link{extractTranscriptSeqs}} for extracting transcript (or CDS)
#' sequences from chromosome sequences.
#' 
#' \item \code{\link{coverageByTranscript}} for computing coverage by
#' transcript (or CDS) of a set of ranges.
#' 
#' \item \code{\link{makeTxDbFromUCSC}}, \code{\link{makeTxDbFromBiomart}}, and
#' \code{\link{makeTxDbFromEnsembl}}, for making a \link{TxDb} object from
#' online resources.
#' 
#' \item \code{\link{makeTxDbFromGRanges}} and \code{\link{makeTxDbFromGFF}}
#' for making a \link{TxDb} object from a \link[GenomicRanges]{GRanges} object,
#' or from a GFF or GTF file.
#' 
#' \item The \link{TxDb} class.  }
#' @keywords manip
#' @examples
#' 
#' library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
#' txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
#' dm3_txlens <- transcriptLengths(txdb)
#' head(dm3_txlens)
#' 
#' dm3_txlens <- transcriptLengths(txdb, with.cds_len=TRUE,
#'                                       with.utr5_len=TRUE,
#'                                       with.utr3_len=TRUE)
#' head(dm3_txlens)
#' 
#' ## When cds_len is 0 (non-coding transcript), utr5_len and utr3_len
#' ## must also be 0:
#' non_coding <- dm3_txlens[dm3_txlens$cds_len == 0, ]
#' stopifnot(all(non_coding[6:8] == 0))
#' 
#' ## When cds_len is not 0 (coding transcript), cds_len + utr5_len +
#' ## utr3_len must be equal to tx_len:
#' coding <- dm3_txlens[dm3_txlens$cds_len != 0, ]
#' stopifnot(all(rowSums(coding[6:8]) == coding[[5]]))
#' 
#' ## A sanity check:
#' stopifnot(identical(dm3_txlens$tx_id, mcols(transcripts(txdb))$tx_id))
#' 
#' @export transcriptLengths
transcriptLengths <- function(txdb, with.cds_len=FALSE,
                                    with.utr5_len=FALSE, with.utr3_len=FALSE,
				    ...)
{
    if (!isTRUEorFALSE(with.cds_len))
        stop("'with.cds_len' must be TRUE or FALSE")
    if (!isTRUEorFALSE(with.utr5_len))
        stop("'with.utr5_len' must be TRUE or FALSE")
    if (!isTRUEorFALSE(with.cds_len))
        stop("'with.utr3_len' must be TRUE or FALSE")
    tx <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"),...)
    ans <- mcols(tx)
    ans$gene_id <- as.character(ans$gene_id)
    tx_id <- as.character(ans$tx_id)  # because match() will want a character

    rg_by_tx <- exonsBy(txdb, by="tx", ...)
    ans$nexon <- .eltNROWS(rg_by_tx, tx_id)
    ans$tx_len <- .sum_width(rg_by_tx, tx_id)
    if (with.cds_len) {
        rg_by_tx <- cdsBy(txdb, by="tx", ...)
        ans$cds_len <- .sum_width(rg_by_tx, tx_id)
    }
    if (with.utr5_len) {
        rg_by_tx <- fiveUTRsByTranscript(txdb, ...)
        ans$utr5_len <- .sum_width(rg_by_tx, tx_id)
    }
    if (with.utr3_len) {
        rg_by_tx <- threeUTRsByTranscript(txdb, ...)
        ans$utr3_len <- .sum_width(rg_by_tx, tx_id)
    }
    as.data.frame(ans)
}

