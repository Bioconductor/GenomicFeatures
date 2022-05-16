#' Map IDs to Genomic Ranges
#'
#'
#' @name mapIdsToRanges
#' @aliases mapIdsToRanges mapIdsToRanges,TxDb-method
#' @docType methods
#' @param x Database to use for mapping
#' @param keys Values to lookup, passed to \code{\link{transcripts}} et. al.
#' @param type Types of feature to return
#' @param columns Additional metadata columns to include in the output
#' @param ... Additional arguments passed to methods
#' @return \code{\link[GenomicRanges]{GRangesList}} corresponding to the keys
#' @section Methods (by class): \itemize{ \item \code{TxDb}: TxDb method }
#' @examples
#'
#' fl <- system.file(package = "GenomicFeatures", "extdata", "sample_ranges.rds")
#' txdb <- makeTxDbFromGRanges(readRDS(fl))
#'
#' keys <- list(tx_name = c("ENST00000371582", "ENST00000371588",
#'     "ENST00000494752", "ENST00000614008", "ENST00000496771"))
#' mapIdsToRanges(txdb, keys = keys, type = "tx")
#'
#' @export
setGeneric("mapIdsToRanges", signature="x",
    function(x, ...) standardGeneric("mapIdsToRanges")
)

setMethod("mapIdsToRanges", "TxDb",
          function(x,
                   keys,
                   type = c("cds", "exon", "tx", "gene"),
                   columns = NULL)
{
    .assert(is.list(keys) && .is.named(keys),
        "'keys' must be a named list")

    .assert(is.null(columns) || is.character(columns),
           "'columns' must be 'NULL' or a character vector")

    type <- match.arg(type)

    fun <- switch(type,
                  cds = cds,
                  exon = exons,
                  tx = transcripts,
                  gene = genes)

    res <- fun(x, keys, columns = unique(c(names(keys), columns)))
    matches <- match(mcols(res)[[names(keys)]], keys[[1]])
    ranges <- rep(res, lengths(matches))

    f <- factor(keys[[1]][unlist(matches, use.names = FALSE)],
                levels = unique(keys[[1]]))
    splitAsList(ranges, f, drop = FALSE)[keys[[1]]]
})

#' Map Genomic Ranges to IDs
#'
#'
#' @name mapRangesToIds
#' @aliases mapRangesToIds mapRangesToIds,TxDb-method
#' @docType methods
#' @param x Database to use for mapping
#' @param ranges range object used to subset
#' @param type of feature to return
#' @param columns additional metadata columns to include in the output.
#' @param ... Additional arguments passed to
#' \code{\link[GenomicRanges:findOverlaps-methods]{findOverlaps}}
#' @return \code{\link[S4Vectors]{DataFrame}} of mcols from the database.
#' @section Methods (by class): \itemize{ \item \code{TxDb}: TxDb method }
#' @examples
#'
#' fl <- system.file(package = "GenomicFeatures", "extdata", "sample_ranges.rds")
#' txdb <- makeTxDbFromGRanges(readRDS(fl))
#'
#' keys <- list(tx_name = c("ENST00000371582", "ENST00000371588",
#'     "ENST00000494752", "ENST00000614008", "ENST00000496771"))
#' res <- mapIdsToRanges(txdb, keys = keys, type = "tx")
#' mapRangesToIds(txdb, res, "tx")
#'
#' @export
setGeneric("mapRangesToIds", signature="x",
    function(x, ...) standardGeneric("mapRangesToIds")
)

setMethod("mapRangesToIds", "TxDb",
          function(x,
                   ranges,
                   type = c("cds", "exon", "tx", "gene"),
                   columns = NULL,
                   ...)
{
    type <- match.arg(type)
    .assert(is(ranges, "Vector"),
        "'ranges' must be a 'Vector'")
    .assert(is.null(columns) || is.character(columns),
           "'columns' must be 'NULL' or a character vector")

    fun <- switch(type,
                  cds = cds,
                  exon = exons,
                  tx = transcripts,
                  gene = genes)

    all <-
        if (is.null(columns)) {
            fun(x)
        } else {
            fun(x, columns = columns)
        }

    hits <- findOverlaps(ranges, all, ...)
    lapply(split(all[subjectHits(hits)], names(ranges)[queryHits(hits)]), mcols)
})

.assert <- function(x, message) {
    if(!x) {
        stop(message, call. = FALSE)
    }
}

.is.named <- function(x) {
    nm <- names(x)
    !is.null(nm) && all(!is.na(nm) & nzchar(nm))
}
