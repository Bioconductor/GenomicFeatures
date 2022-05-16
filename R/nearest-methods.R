### =========================================================================
### nearest (and related) methods
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### distance
###

#' Finding the nearest genomic range neighbor in a TxDb
#'
#' The \code{distance} methods for TxDb objects and subclasses.
#'
#' \itemize{ \itemdistance: Returns the distance for each range in \code{x} to
#' the range extracted from the \link{TxDb} object \code{y}. Values in
#' \code{id} are matched to one of \sQuote{gene_id}, \sQuote{tx_id},
#' \sQuote{exon_id} or \sQuote{cds_id} identifiers in the \link{TxDb} and the
#' corresponding ranges are extracted. The \code{type} argument specifies which
#' identifier is represented in \code{id}. The extracted ranges are used in the
#' distance calculation with the ranges in \code{x}.
#'
#' The method returns \code{NA} values when the genomic region defined by
#' \code{id} cannot be collapsed into a single range (e.g., when a gene spans
#' multiple chromosomes) or if the \code{id} is not found in \code{y}.
#'
#' The behavior of \code{distance} with respect to zero-width ranges has
#' changed in Bioconductor 2.12. See the man page \code{?distance} in IRanges
#' for details.
#'
#' }
#'
#' @aliases nearest-methods distance,GenomicRanges,TxDb-method
#' @param x The query \link{GenomicRanges} instance.
#' @param y For \code{distance}, a \link{TxDb} instance. The \code{id} is used
#' to extract ranges from the \link{TxDb} which are then used to compute the
#' distance from \code{x}.
#' @param id A \code{character} vector the same length as \code{x}.  The
#' \code{id} must be identifiers in the \link{TxDb} object.  \code{type}
#' indicates what type of identifier \code{id} is.
#' @param type A \code{character(1)} describing the \code{id}.  Must be one of
#' \sQuote{gene}, \sQuote{tx}, \sQuote{exon} or \sQuote{cds}.
#' @param ignore.strand A \code{logical} indicating if the strand of the ranges
#' should be ignored. When \code{TRUE}, strand is set to \code{'+'}.
#' @param ... Additional arguments for methods.
#' @return For \code{distance}, an integer vector of distances between the
#' ranges in \code{x} and \code{y}.
#' @author Valerie Obenchain <vobencha@@fhcrc.org>
#' @seealso \itemize{ \item \link[IRanges]{nearest-methods} man page in
#' IRanges.  \item \link[GenomicRanges]{nearest-methods} man page in
#' GenomicRanges.  }
#' @keywords utilities
#' @examples
#'
#'   ## -----------------------------------------------------------
#'   ## distance()
#'   ## -----------------------------------------------------------
#'
#'   library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
#'   txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
#'   gr <- GRanges(c("chr2L", "chr2R"),
#'                 IRanges(c(100000, 200000),  width=100))
#'   distance(gr, txdb, id=c("FBgn0259717", "FBgn0261501"), type="gene")
#'   distance(gr, txdb, id=c("10000", "23000"), type="cds")
#'
#'   ## The id's must be in the appropriate order with respect to 'x'.
#'   distance(gr, txdb, id=c("4", "4097"), type="tx")
#'
#'   ## 'id' "4" is on chr2L and "4097" is on chr2R.
#'   transcripts(txdb, filter=list(tx_id=c("4", "4097")))
#'
#'   ## If we reverse the 'id' the chromosomes are incompatable with gr.
#'   distance(gr, txdb, id=c("4097", "4"), type="tx")
#'
#'   ## distance() compares each 'x' to the corresponding 'y'.
#'   ## If an 'id' is not found in the TxDb 'y' will not
#'   ## be the same lenth as 'x' and an error is thrown.
#'   \dontrun{
#'   distance(gr, txdb, id=c("FBgn0000008", "INVALID"), type="gene") ## will fail
#'   }
#'
#' @export
setMethod("distance", c("GenomicRanges", "TxDb"),
    function(x, y, ignore.strand=FALSE, ..., id,
             type=c("gene", "tx", "exon", "cds"))
    {
        if (!identical(length(x), length(id)))
            stop("length(id) must equal length(x)")
        if (!is.character(id))
            stop("'id' must be a character")

        if (type == "gene") {
            .extractByGeneID(x, y, ignore.strand, id)
        } else {
            rng <- switch(type,
                          tx=transcripts(y, "tx_id", filter=list(tx_id=id)),
                          exon=exons(y, "exon_id", filter=list(exon_id=id)),
                          cds=cds(y, "cds_id", filter=list(cds_id=id)))
            f <- factor(mcols(rng)[,])
            missing <- !id %in% levels(f)
            if (any(missing))
                  warning(paste0("id(s): '", paste(unique(id[missing]),
                                 sep=","), "' were not found in 'y'"))
            ## rep out ranges according to 'id'
            rng <- rng[match(id[!missing], levels(f))]
            ans <- rep(NA_integer_, length(x))
            ans[!missing] <- distance(x[!missing], rng,
                                      ignore.strand=ignore.strand)
            stopifnot(length(ans) == length(id))
            ans
        }
    }
)

.extractByGeneID <- function(x, y, ignore.strand, id)
{
    tx <- transcriptsBy(y, "gene")
    missing <- !id %in% names(tx)
    if (any(missing))
          warning(paste0("id(s): '", paste(unique(id[missing]), sep=","),
                  "' were not found in 'y'"))

    group <- range(tx[names(tx) %in% id], ignore.strand=ignore.strand)
    multiRange <- lengths(group) > 1L
    if (any(multiRange)) {
        warning(paste0("id(s): '", paste(unique(names(multiRange)[multiRange]),
                       sep=','),
                       "' could not be collapsed to a single gene region"))
        group <- group[!multiRange]
    }

    valid <- (!id %in% names(multiRange)[multiRange]) & !missing
    ## rep out ranges according to 'id'
    rng <- unlist(group, use.names=FALSE)
    rng <- rng[match(id[valid], names(group))]
    ans <- rep(NA_integer_, length(x))
    ans[valid] <- distance(x[valid], rng, ignore.strand=ignore.strand)
    stopifnot(length(ans) == length(id))
    ans
}
