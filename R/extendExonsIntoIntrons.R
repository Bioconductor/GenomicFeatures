### =========================================================================
### Extend exons by a given number of bases into their adjacent introns
### -------------------------------------------------------------------------
###



#' Extend exons by a given number of bases into their adjacent introns
#' 
#' \code{extendExonsIntoIntrons} extends the supplied exons by a given number
#' of bases into their adjacent introns.
#' 
#' 
#' @param ex_by_tx A \link[GenomicRange]{GRangesList} object containing exons
#' grouped by transcript. This must be an object as returned by
#' \code{\link{exonsBy}(txdb, by="tx")}, that is: \itemize{ \item each list
#' element in \code{ex_by_tx} must be a \link[GenomicRange]{GRanges} object
#' representing the exons of a given transcript; \item the exons in each list
#' element must be ordered by ascending rank with respect to their transcript.
#' }
#' @param extent Size of the extent in number of bases. 2 by default.
#' 
#' The first exon in a transcript will be extended by that amount on its 3'
#' side only. The last exon in a transcript will be extended by that amount on
#' its 5' side only. All other exons (i.e. intermediate exons) will be extended
#' by that amount on \emph{each} side.
#' 
#' Note that exons that belong to a single-exon transcript don't get extended.
#' 
#' The default value of 2 corresponds to inclusion of the donor/acceptor
#' intronic regions (typically GT/AG).
#' @return A copy of \link[GenomicRange]{GRangesList} object \code{ex_by_tx}
#' where the original exon ranges have been extended.
#' 
#' Names and metadata columns on \code{ex_by_tx} are propagated to the result.
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
#' \item The \link{TxDb} class.  }
#' @keywords manip
#' @examples
#' 
#' ## With toy transcripts:
#' ex_by_tx <- GRangesList(
#'     TX1="chr1:10-20:+",
#'     TX2=c("chr1:10-20:+", "chr1:50-75:+"),
#'     TX3=c("chr1:10-20:+", "chr1:50-75:+", "chr1:100-120:+"),
#'     TX4="chr1:10-20:-",
#'     TX5=c("chr1:10-20:-", "chr1:50-75:-"),
#'     TX6=c("chr1:10-20:-", "chr1:50-75:-", "chr1:100-120:-")
#' )
#' 
#' extended <- extendExonsIntoIntrons(ex_by_tx, extent=2)
#' extended[1:3]
#' extended[4:6]
#' 
#' ## With real-world transcripts:
#' library(TxDb.Celegans.UCSC.ce11.ensGene)
#' txdb <- TxDb.Celegans.UCSC.ce11.ensGene
#' ex_by_tx <- exonsBy(txdb, by="tx")
#' ex_by_tx
#' 
#' extendExonsIntoIntrons(ex_by_tx, extent=2)
#' 
#' ## Sanity check:
#' stopifnot(identical(extendExonsIntoIntrons(ex_by_tx, extent=0), ex_by_tx))
#' 
#' @export extendExonsIntoIntrons
extendExonsIntoIntrons <- function(ex_by_tx, extent=2)
{
    if (!is(ex_by_tx, "GRangesList"))
        stop(wmsg("'ex_by_tx' must be a GRangesList object"))
    if (!isSingleNumber(extent))
        stop(wmsg("'extent' must be a single number"))
    if (!is.integer(extent))
        extent <- as.integer(extent)

    resize_idx <- which(lengths(ex_by_tx) >= 2L)
    ex_to_resize <- ex_by_tx[resize_idx]

    ## Resize first exons.
    first_ex <- heads(ex_to_resize, n=1L)
    unlisted <- unlist(first_ex, use.names=FALSE)
    unlisted <- resize(unlisted, width(unlisted) + extent,
                       fix="start", use.names=FALSE)
    first_ex <- relist(unlisted, first_ex)

    ## Resize last exons.
    last_ex <- tails(ex_to_resize, n=1L)
    unlisted <- unlist(last_ex, use.names=FALSE)
    unlisted <- resize(unlisted, width(unlisted) + extent,
                       fix="end", use.names=FALSE)
    last_ex <- relist(unlisted, last_ex)

    ## Resize intermediate exons.
    mid_ex <- tails(heads(ex_to_resize, n=-1L), n=-1L)
    unlisted <- unlist(mid_ex, use.names=FALSE)
    unlisted <- resize(unlisted, width=width(unlisted) + 2L*extent,
                       fix="center", use.names=FALSE)
    mid_ex <- relist(unlisted, mid_ex)

    ## Put exons back together.
    ex_by_tx[resize_idx] <- pc(first_ex, mid_ex, last_ex)
    ex_by_tx
}

