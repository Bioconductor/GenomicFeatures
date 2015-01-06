### =========================================================================
### mapToTranscript() and pmapToTranscript() methods
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generics 
###

setGeneric("mapToTranscript", signature=c("x", "transcript"),
    function(x, transcript, reverse=FALSE, ...) 
        standardGeneric("mapToTranscript")
)

setGeneric("pmapToTranscript", signature=c("x", "transcript"),
    function(x, transcript, reverse=FALSE, ...) 
        standardGeneric("pmapToTranscript")
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
### mapToTranscript()
###

### reverse = FALSE
.mapToTranscript <- function(x, transcript, hits, ignore.strand) 
{
    flat <- unlist(transcript)
    if (length(hits)) {
        xHits <- queryHits(hits)
        transcriptHits <- subjectHits(hits)
        xrange <- ranges(x)[xHits]
        bounds <- ranges(flat)[transcriptHits]

        ## location wrt to start of individual list elements
        if (ignore.strand) {
            xrange <- shift(xrange, - start(bounds) + 1L)
        } else {
            neg <- as.vector(strand(flat)[transcriptHits] == "-")
            negstart <- end(bounds)[neg] - end(xrange)[neg] + 1L
            xrange[neg] <- IRanges(negstart, width=width(xrange)[neg])
            xrange[!neg] <- shift(xrange, - start(bounds) + 1L)
        }
        ## location wrt start of concatenated list elements
        if (length(flat) > length(transcript)) {
            shifted <- .listCumsumShifted(width(transcript))
            xrange <- shift(xrange, shifted[transcriptHits])
        }
        GRanges(seqnames(flat)[transcriptHits], xrange, 
                strand(flat)[transcriptHits],
                DataFrame(xHits, 
                          transcriptHits=togroup(transcript)[transcriptHits]))
    } else {
        ans <- GRanges()
        mcols(ans) <- DataFrame(xHits=integer(), transcriptHits=integer())
        ans
    }
}

### reverse = TRUE 
.rmapToTranscript <- function(x, transcript, hits, ignore.strand)
{
    if (length(hits)) {
        xHits <- queryHits(hits)
        transcriptHits <- subjectHits(hits)
        ## strand mismatch
        if (!ignore.strand) {
            txstrand <- unlist(runValue(strand(transcript)), use.names=FALSE)
            mismatch <- (as.numeric(strand(x)[xHits]) + 
                         as.numeric(txstrand[transcriptHits])) == 3L 
            strand <- as.character(txstrand)[transcriptHits]
        } else {
            mismatch <- logical(length(xHits))
            strand <- rep("+", length(xHits))
        }
        xStart <- as.list(start(x)[xHits])
        xEnd <- as.list(end(x)[xHits])
        transcriptStart <- as.list(start(transcript)[transcriptHits])
        transcriptEnd <- as.list(end(transcript)[transcriptHits])
        s <- unlist(transcriptLocs2refLocs(xStart, transcriptStart, 
                                           transcriptEnd, strand,
                                           FALSE, FALSE), use.names=FALSE)
        e <- unlist(transcriptLocs2refLocs(xEnd, transcriptStart, 
                                           transcriptEnd, strand,
                                           FALSE, FALSE), use.names=FALSE)

        ## non-hits are zero width with seqname "unmmapped"
        seqname <- as.character(runValue(seqnames(transcript)[transcriptHits]))
        if (any(skip <- is.na(s) | is.na(e) | mismatch)) {
            s[skip] <- 1L 
            e[skip] <- 0L
            seqname[skip] <- "unmapped"
        }
        GRanges(Rle(seqname), IRanges(s, e, names=names(x)[xHits]), 
                strand=strand, DataFrame(xHits, transcriptHits))
    } else {
        ans <- GRanges()
        mcols(ans) <- DataFrame(xHits=integer(), transcriptHits=integer())
        ans
    }
}

setMethod("mapToTranscript", c("GenomicRanges", "GenomicRanges"), 
    function(x, transcript, reverse=FALSE, ignore.strand=TRUE, ...)
    {
        grl <- relist(transcript, PartitioningByEnd(seq_along(transcript), 
                      names=names(transcript)))
        if (reverse)
            mapToTranscript(x, grl, TRUE, ignore.strand, ...)
        else 
            mapToTranscript(x, grl, FALSE, ignore.strand, ...)
    }
)

setMethod("mapToTranscript", c("GenomicRanges", "GRangesList"), 
    function(x, transcript, reverse=FALSE, ignore.strand=TRUE, ...) 
    {
        if (length(x) && length(transcript)) {
            if (!ignore.strand)
                if (!all(elementLengths(runLength(strand(transcript))) == 1))
                    stop(paste0("when ignore.strand=TRUE all inner list ",
                                "elements of 'transcript' must be the ",
                                "same strand"))

            ## order within list elements by strand
            transcript <- 
                .orderElementsByTranscription(transcript, ignore.strand)
            if (reverse) {
                ## name matching determines pairs
                if (is.null(xNames <- names(x)) || 
                    is.null(transcriptNames <- names(transcript)))
                    stop ("both 'x' and 'transcript' must have names")
                match0 <- match(transcriptNames, transcriptNames)
                match1 <- match(xNames, transcriptNames)
                group0 <- splitAsList(seq_along(transcriptNames), match0)
                group1 <- group0[match(na.omit(match1), names(group0))]
                xHits <- rep(which(!is.na(match1)), elementLengths(group1))
                transcriptHits <- unlist(group1, use.names=FALSE)
                if (!length(xHits <- na.omit(xHits)))
                    stop ("none of 'names(x)' are in 'names(transcript)'")

                hits <- Hits(xHits, transcriptHits, length(x), 
                             length(transcript))
                map <- .rmapToTranscript(x, transcript, hits, ignore.strand) 
                ## remove zero-width ranges 
                if (any(hasWidth <- width(map) != 0L))
                    map <- map[hasWidth]
                map
            } else {
                ## findOverlaps determines pairs
                hits <- findOverlaps(x, unlist(transcript, use.names=FALSE), 
                                     type="within", ignore.strand=ignore.strand)
                .mapToTranscript(x, transcript, hits, ignore.strand)

            }
        } else {
            ans <- GRanges()
            mcols(ans) <- DataFrame(xHits=integer(), transcriptHits=integer())
            ans
        }
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pmapToTranscript()
###

setMethod("pmapToTranscript", c("Ranges", "GenomicRanges"),
    function(x, transcript, reverse=FALSE, ...) 
    { 
        gr <- GRanges(seqnames(transcript), x)
        if (reverse)
            ranges(pmapToTranscript(gr, transcript, TRUE, TRUE))
        else 
            ranges(pmapToTranscript(gr, transcript, FALSE, TRUE))
    }
)

setMethod("pmapToTranscript", c("GenomicRanges", "GenomicRanges"), 
    function(x, transcript, reverse=FALSE,  ignore.strand=TRUE, ...) 
    {
        grl <- splitAsList(transcript, seq_along(transcript))
        if (reverse)
            pmapToTranscript(x, grl, TRUE, ignore.strand, ...)
        else 
            pmapToTranscript(x, grl, FALSE, ignore.strand, ...)
    }
)

setMethod("pmapToTranscript", c("GenomicRanges", "GRangesList"), 
    function(x, transcript, ignore.strand=TRUE, ...) 
    {
        if (length(x) && length(transcript)) {
            if (length(x) != length(transcript))
                stop("'x' and 'transcript' must have the same length")
            if (!ignore.strand)
                if (!all(elementLengths(runLength(strand(transcript))) == 1))
                    stop(paste0("when ignore.strand=TRUE all inner list ",
                                "elements of 'transcript' must have the ",
                                "same strand"))
            ## order within list elements
            transcript <- 
                .orderElementsByTranscription(transcript, ignore.strand)

            if (reverse) {
                ## i-th element matching determines pairs
                hits <- Hits(seq_along(x), seq_along(x), length(x), length(x))
                .rmapToTranscript(x, transcript, hits, ignore.strand)
            } else {
                ## i-th element matching determines pairs
                hits <- findOverlaps(x, unlist(transcript, use.names=FALSE), 
                                     type="within", ignore.strand=ignore.strand)
                ith <- 
                    queryHits(hits) == togroup(transcript)[subjectHits(hits)]
                map <- .mapToTranscript(x, transcript, hits[ith], ignore.strand)

                ## non-hits are zero width with seqname "unmmapped"
                if (length(x) != length(map)) {
                    starts <- rep(1L, length(x))
                    ends <- rep(0L, length(x))
                    seqname <- Rle("unmapped", length(x))
                    strands <- Rle("*", length(x))
                    xHits <- mcols(map)$xHits
                    ends[xHits] <- end(map)
                    starts[xHits] <- start(map)
                    seqname[xHits] <- seqnames(map)
                    strands[xHits] <- strand(map)
                    GRanges(seqname, IRanges(starts, ends, names=names(x)), 
                            strands)
                } else {
                    mcols(map) <- NULL
                    map
                }
            }
        } else {
            GRanges(rep("unmapped", length(x)), 
                    IRanges(rep(1, length(x)), width=0))
        }
    }
)

