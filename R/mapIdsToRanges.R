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
