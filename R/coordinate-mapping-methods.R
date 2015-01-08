### =========================================================================
### mapToTranscripts() and pmapToTranscripts() methods
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generics 
###

setGeneric("mapToTranscripts", signature=c("x", "transcripts"),
    function(x, transcripts, reverse=FALSE, ...) 
        standardGeneric("mapToTranscripts")
)

setGeneric("pmapToTranscripts", signature=c("x", "transcripts"),
    function(x, transcripts, reverse=FALSE, ...) 
        standardGeneric("pmapToTranscripts")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers
###

### 'x' is a GRangesList
### Returns a GRangesList with sorted elements. This method differs from 
### sort() in that "-" strand elements are returned highest value to lowest.
.orderElementsByTranscription <- function(x, ignore.strand) {
    original <- unlist(sapply(elementLengths(x), function(xx) 1:xx), 
                       use.names=FALSE)
    ## order by position
    gr <- unlist(x, use.names = FALSE)
    idx <- order(togroup(x), start(gr))
    gr <- gr[idx]
    part <- PartitioningByWidth(x)
    ## handle zero-width ranges
    pstart <- start(part)[width(part) != 0L]
    pend <- end(part)[width(part) != 0L]

    if (ignore.strand) {
        ord <- S4Vectors:::mseq(pstart, pend)
    } else {
        neg <- strand(gr)[pstart] == "-"
        ord <- S4Vectors:::mseq(ifelse(neg, pend, pstart),
                                ifelse(neg, pstart, pend))
    }
    res <- relist(gr[ord], x)
    res@unlistData$unordered <- original[idx[ord]] 
    res
}

### 'x' is an IntegerList or NumericList
### Returns a numeric vector of cumulative sums within list elements.
.listCumsumShifted <- function(x) {
    cs <- unlist(cumsum(x), use.names=FALSE)
    shifted <- c(0L, head(cs, -1))
    shifted[start(PartitioningByWidth(elementLengths(x)))] <- 0L
    shifted
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToTranscripts()
###

### reverse = FALSE
.mapToTranscripts <- function(x, transcripts, hits, ignore.strand) 
{
    flat <- unlist(transcripts)
    if (length(hits)) {
        xHits <- queryHits(hits)
        txHits <- subjectHits(hits)
        xrange <- ranges(x)[xHits]
        bounds <- ranges(flat)[txHits]

        ## location wrt to start of individual list elements
        if (ignore.strand) {
            xrange <- shift(xrange, - start(bounds) + 1L)
        } else {
            neg <- as.vector(strand(flat)[txHits] == "-")
            negstart <- end(bounds)[neg] - end(xrange)[neg] + 1L
            xrange[neg] <- IRanges(negstart, width=width(xrange)[neg])
            xrange[!neg] <- shift(xrange, - start(bounds) + 1L)
        }
        ## location wrt start of concatenated list elements
        if (length(flat) > length(transcripts)) {
            shifted <- .listCumsumShifted(width(transcripts))
            xrange <- shift(xrange, shifted[txHits])
        }
        ## seqnames come from 'transcripts'
        GRanges(rep(names(flat), elementLengths(flat))[txHits],
                xrange, strand(flat)[txHits],
                DataFrame(xHits, 
                          transcriptsHits=togroup(transcripts)[txHits])) 
    } else {
        ans <- GRanges()
        mcols(ans) <- DataFrame(xHits=integer(), transcriptsHits=integer())
        ans
    }
}

### reverse = TRUE 
.rmapToTranscripts <- function(x, transcripts, hits, ignore.strand)
{
    if (length(hits)) {
        xHits <- queryHits(hits)
        txHits <- subjectHits(hits)
        ## strand mismatch
        if (!ignore.strand) {
            txstrand <- unlist(runValue(strand(transcripts)), use.names=FALSE)
            mismatch <- (as.numeric(strand(x)[xHits]) + 
                         as.numeric(txstrand[txHits])) == 3L 
            strand <- as.character(txstrand)[txHits]
        } else {
            mismatch <- logical(length(xHits))
            strand <- rep("+", length(xHits))
        }
        xStart <- as.list(start(x)[xHits])
        xEnd <- as.list(end(x)[xHits])
        txStart <- as.list(start(transcripts)[txHits])
        txEnd <- as.list(end(transcripts)[txHits])
        s <- unlist(transcriptLocs2refLocs(xStart, txStart, txEnd, strand,
                                           FALSE, FALSE), use.names=FALSE)
        e <- unlist(transcriptLocs2refLocs(xEnd, txStart, txEnd, strand,
                                           FALSE, FALSE), use.names=FALSE)

        ## non-hits are zero width with seqname "unmmapped"
        seqname <- as.character(runValue(seqnames(transcripts)[txHits]))
        if (any(skip <- is.na(s) | is.na(e) | mismatch)) {
            s[skip] <- 1L 
            e[skip] <- 0L
            seqname[skip] <- "unmapped"
        }
        GRanges(Rle(seqname), IRanges(s, e, names=names(x)[xHits]), 
                strand=strand, DataFrame(xHits, transcriptsHits=txHits))
    } else {
        ans <- GRanges()
        mcols(ans) <- DataFrame(xHits=integer(), transcriptsHits=integer())
        ans
    }
}

setMethod("mapToTranscripts", c("GenomicRanges", "GenomicRanges"), 
    function(x, transcripts, reverse=FALSE, ignore.strand=TRUE, ...)
    {
        grl <- relist(transcripts, PartitioningByEnd(seq_along(transcripts), 
                      names=names(transcripts)))
        if (reverse)
            mapToTranscripts(x, grl, TRUE, ignore.strand, ...)
        else 
            mapToTranscripts(x, grl, FALSE, ignore.strand, ...)
    }
)

setMethod("mapToTranscripts", c("GenomicRanges", "GRangesList"), 
    function(x, transcripts, reverse=FALSE, ignore.strand=TRUE, ...) 
    {
        if (length(x) && length(transcripts)) {
            if (!ignore.strand)
                if (!all(elementLengths(runLength(strand(transcripts))) == 1))
                    stop(paste0("when ignore.strand=TRUE all inner list ",
                                "elements of 'transcripts' must be the ",
                                "same strand"))

            ## order within list elements by strand
            transcripts <- 
                .orderElementsByTranscription(transcripts, ignore.strand)
            if (reverse) {
                ## name matching determines pairs
                if (is.null(xNames <- names(x)) || 
                    is.null(transcriptsNames <- names(transcripts)))
                    stop ("both 'x' and 'transcripts' must have names")
                match0 <- match(transcriptsNames, transcriptsNames)
                match1 <- match(xNames, transcriptsNames)
                group0 <- splitAsList(seq_along(transcriptsNames), match0)
                group1 <- group0[match(na.omit(match1), names(group0))]
                xHits <- rep(which(!is.na(match1)), elementLengths(group1))
                txHits <- unlist(group1, use.names=FALSE)
                if (!length(xHits <- na.omit(xHits)))
                    stop ("none of 'names(x)' are in 'names(transcripts)'")

                hits <- Hits(xHits, txHits, length(x), length(transcripts))
                map <- .rmapToTranscripts(x, transcripts, 
                                            hits, ignore.strand) 
                ## remove zero-width ranges 
                if (any(hasWidth <- width(map) != 0L))
                    map <- map[hasWidth]
                map
            } else {
                if (is.null(names(transcripts)))
                    stop ("'transcripts' must have names")
                ## findOverlaps determines pairs
                hits <- findOverlaps(x, unlist(transcripts, use.names=FALSE), 
                                     type="within", ignore.strand=ignore.strand)
                .mapToTranscripts(x, transcripts, hits, ignore.strand)

            }
        } else {
            ans <- GRanges()
            mcols(ans) <- DataFrame(xHits=integer(), 
                                    transcriptsHits=integer())
            ans
        }
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pmapToTranscripts()
###

setMethod("pmapToTranscripts", c("Ranges", "GenomicRanges"),
    function(x, transcripts, reverse=FALSE, ...) 
    { 
        gr <- GRanges(seqnames(transcripts), x)
        if (reverse)
            ranges(pmapToTranscripts(gr, transcripts, TRUE, TRUE))
        else 
            ranges(pmapToTranscripts(gr, transcripts, FALSE, TRUE))
    }
)

setMethod("pmapToTranscripts", c("GenomicRanges", "GenomicRanges"), 
    function(x, transcripts, reverse=FALSE,  ignore.strand=TRUE, ...) 
    {
        grl <- splitAsList(transcripts, seq_along(transcripts))
        if (reverse)
            pmapToTranscripts(x, grl, TRUE, ignore.strand, ...)
        else 
            pmapToTranscripts(x, grl, FALSE, ignore.strand, ...)
    }
)

setMethod("pmapToTranscripts", c("GenomicRanges", "GRangesList"), 
    function(x, transcripts, reverse=FALSE, ignore.strand=TRUE, ...) 
    {
        if (length(x) && length(transcripts)) {
            if (length(x) != length(transcripts))
                stop("'x' and 'transcripts' must have the same length")
            if (!ignore.strand)
                if (!all(elementLengths(runLength(strand(transcripts))) == 1))
                    stop(paste0("when ignore.strand=TRUE all inner list ",
                                "elements of 'transcripts' must have the ",
                                "same strand"))
            ## order within list elements
            transcripts <- 
                .orderElementsByTranscription(transcripts, ignore.strand)

            if (reverse) {
                ## i-th element matching determines pairs
                hits <- Hits(seq_along(x), seq_along(x), length(x), length(x))
                .rmapToTranscripts(x, transcripts, hits, ignore.strand)
            } else {
                if (is.null(names(transcripts)))
                    stop ("'transcripts' must have names")
                ## i-th element matching determines pairs
                hits <- findOverlaps(x, unlist(transcripts, use.names=FALSE), 
                                     type="within", ignore.strand=ignore.strand)
                ith <- 
                    queryHits(hits) == togroup(transcripts)[subjectHits(hits)]
                map <- .mapToTranscripts(x, transcripts, 
                                           hits[ith], ignore.strand)

                ## non-hits are zero width with seqname "unmmapped"
                if (length(x) != length(map)) {
                    s <- rep(1L, length(x))
                    e <- rep(0L, length(x))
                    seqname <- Rle("unmapped", length(x))
                    strands <- Rle("*", length(x))
                    xHits <- mcols(map)$xHits
                    e[xHits] <- end(map)
                    s[xHits] <- start(map)
                    seqname[xHits] <- seqnames(map)
                    strands[xHits] <- strand(map)
                    GRanges(seqname, IRanges(s, e, names=names(x)), 
                            strands)
                } else {
                    mcols(map) <- NULL
                    map
                }
            }
        } else {
            GRanges(rep("unmapped", length(x)), 
                    IRanges(rep(1, length(x)), width=0, names=names(x)))
        }
    }
)

