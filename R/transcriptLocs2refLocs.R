### =========================================================================
### transcriptLocs2refLocs()
### -------------------------------------------------------------------------


.normargExonStartsOrEnds <- function(exonStarts, argname)
{
    if (is.list(exonStarts))
        return(exonStarts)
    if (is(exonStarts, "IntegerList"))
        return(as.list(exonStarts))
    if (is.character(exonStarts))
        return(strsplitAsListOfIntegerVectors(exonStarts))
    stop("'", argname, "' must be a list of integer vectors, ",
         "an IntegerList object,\n  or a character vector where ",
         "each element is a comma-separated list of\n  integers")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### transcriptLocs2refLocs()
###

transcriptLocs2refLocs <- function(tlocs, exonStarts=list(), exonEnds=list(),
                                   strand=character(0),
                                   decreasing.rank.on.minus.strand=FALSE,
                                   error.if.out.of.bounds=TRUE)
{
    if (!is.list(tlocs)) {
        if (!is(tlocs, "IntegerList"))
            stop("'tlocs' must be a list of integer vectors ",
                 "or an IntegerList object")
        tlocs <- as.list(tlocs)
    }
    if (is(exonStarts, "RangesList")) {
        if (!identical(exonEnds, list()))
            stop("'exonEnds' cannot be specified ",
                 "when 'exonStarts' is a RangesList object")
        exonEnds <- end(exonStarts)
        exonStarts <- start(exonStarts)
    }
    exonStarts <- .normargExonStartsOrEnds(exonStarts, "exonStarts")
    exonEnds <- .normargExonStartsOrEnds(exonEnds, "exonEnds")
    if (is.factor(strand))
        strand <- as.vector(strand)
    if (!is.character(strand))
        stop("'strand' must be a character vector")
    if (length(tlocs) != length(strand)
     || length(exonStarts) != length(strand)
     || length(exonEnds) != length(strand))
        stop("'tlocs', 'exonStarts', 'exonEnds' and 'strand' ",
             "must have the same length")
    if (!isTRUEorFALSE(decreasing.rank.on.minus.strand))
        stop("'decreasing.rank.on.minus.strand' must be TRUE or FALSE")
    GenomicRanges:::unsafe.transcriptLocs2refLocs(tlocs,
                            exonStarts, exonEnds, strand,
                            decreasing.rank.on.minus.strand,
                            error.if.out.of.bounds)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### transcriptWidths()
###

transcriptWidths <- function(exonStarts=list(), exonEnds=list())
{
    if (is(exonStarts, "RangesList")) {
        if (!identical(exonEnds, list()))
            stop("'exonEnds' cannot be specified ",
                 "when 'exonStarts' is a RangesList object")
        exonEnds <- end(exonStarts)
        exonStarts <- start(exonStarts)
    }
    exonStarts <- .normargExonStartsOrEnds(exonStarts, "exonStarts")
    exonEnds <- .normargExonStartsOrEnds(exonEnds, "exonEnds")
    if (length(exonStarts) != length(exonEnds))
        stop("'exonStarts', 'exonEnds' must have the same length")
    GenomicRanges:::unsafe.transcriptWidths(exonStarts, exonEnds)
}

