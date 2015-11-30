setGeneric("mapIdsToRanges", signature="db",
    function(db, ...) standardGeneric("mapIdsToRanges")
)

setMethod("mapIdsToRanges", "TxDb",
          function(db,
                   vals,
                   type = c("cds", "exon", "tx", "gene"),
                   columns = NULL,
                   ...)
{
    assert(is.list(vals) && is.named(vals),
        "'vals' must be a named list")

    assert(is.null(columns) || is.character(columns),
           "'columns' must be 'NULL' or a character vector")

    type <- match.arg(type)

    fun <- switch(type,
                  cds = cds,
                  exon = exons,
                  tx = transcripts,
                  gene = genes)

    res <- fun(db, vals, columns = unique(c(names(vals), columns)))
    matches <- match(mcols(res)[[names(vals)]], vals[[1]])
    ranges <- rep(res, lengths(matches))

    splitAsList(ranges, factor(vals[[1]][unlist(matches, use.names = FALSE)], levels = vals[[1]]), drop = FALSE)
})

setGeneric("mapRangesToIds", signature="db",
    function(db, ...) standardGeneric("mapRangesToIds")
)

setMethod("mapRangesToIds", "TxDb",
          function(db,
                   ranges,
                   type = c("cds", "exon", "tx", "gene"),
                   columns = NULL,
                   ...)
{
    type <- match.arg(type)
    assert(is(ranges, "Vector"),
        "'ranges' must be a 'Vector'")
    assert(is.null(columns) || is.character(columns),
           "'columns' must be 'NULL' or a character vector")

    fun <- switch(type,
                  cds = cds,
                  exon = exons,
                  tx = transcripts,
                  gene = genes)

    all <-
        if (is.null(columns)) {
            fun(db)
        } else {
            fun(db, columns = columns)
        }

    hits <- findOverlaps(ranges, all, ...)
    lapply(split(all[subjectHits(hits)], names(ranges)[queryHits(hits)]), mcols)
})

assert <- function(x, message) {
    if(!x) {
        stop(message, call. = FALSE)
    }
}

is.named <- function(x) {
    nm <- names(x)
    !is.null(nm) && all(!is.na(nm) & nm != "")
}
