### =========================================================================
### getPromoterSeq() and getTerminatorSeq()
### -------------------------------------------------------------------------

### Original author: Paul Shannon

### NOTE (H. Pagès, Jan 22, 2024): Interface is inconsistent with
### extractTranscriptSeqs() or extractUpstreamSeqs().
### TODO: Implement extractPromoterSeqs() and extractTerminatorSeqs() and
### model them after extractTranscriptSeqs() or extractUpstreamSeqs().
### Then deprecate getPromoterSeq() and getTerminatorSeq() in favor of
### extractPromoterSeqs() and extractTerminatorSeqs().

setGeneric("getPromoterSeq", signature="query",
    function(query, subject, upstream=2000, downstream=200)
        standardGeneric("getPromoterSeq"))

setGeneric("getTerminatorSeq", signature="query",
    function(query, subject, upstream=2000, downstream=200)
        standardGeneric("getTerminatorSeq"))

.GRanges_getPromoterSeq <- function(query, subject, FUN, upstream, downstream)
{
    stopifnot(is(query, "GRanges"))
    promoter.granges <- FUN(query, upstream, downstream)
    result <- getSeq(subject, promoter.granges)
    md <- mcols(query)
    geneIDs <- names(query)   # often NULL
    if (is.null(geneIDs))
      geneIDs <- rep(NA_character_, length(query))
    md$geneID <- geneIDs
    mcols(result) <- md
    result
}

setMethod("getPromoterSeq", "GRanges",
    function(query, subject, upstream=2000, downstream=200)
        .GRanges_getPromoterSeq(query, subject, promoters,
                                upstream, downstream)
)

setMethod("getTerminatorSeq", "GRanges",
    function(query, subject, upstream=2000, downstream=200)
        .GRanges_getPromoterSeq(query, subject, terminators,
                                upstream, downstream)
)

.GRangesList_getPromoterSeq <-
    function(query, subject, FUN, upstream, downstream)
{
    stopifnot(is(query, "GRangesList"))
    unlisted_query <- unlist(query, use.names=FALSE)  # GRanges object
    promoter.granges <- FUN(unlisted_query, upstream, downstream)
    result <- getSeq(subject, promoter.granges)
    md <- mcols(unlisted_query)
    geneIDs <- names(query)
    geneID.counts <- elementNROWS(query)
    geneIDs <- rep(geneIDs, geneID.counts)  # H. Pagès: what if geneIDs is NULL?
    md$geneID <- geneIDs
    mcols(result) <- md
    relist(result, query)
}

setMethod("getPromoterSeq", "GRangesList",
    function(query, subject, upstream=2000, downstream=200)
        .GRangesList_getPromoterSeq(query, subject, promoters,
                                    upstream, downstream)
)

setMethod("getTerminatorSeq", "GRangesList",
    function(query, subject, upstream=2000, downstream=200)
        .GRangesList_getPromoterSeq(query, subject, terminators,
                                    upstream, downstream)
)

