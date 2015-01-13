### =========================================================================
### mapToTranscripts() and pmapToTranscripts() methods
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generics 
###

setGeneric("mapToTranscripts", signature=c("x", "transcripts"),
    function(x, transcripts, ...) 
        standardGeneric("mapToTranscripts")
)

setGeneric("pmapToTranscripts", signature=c("x", "transcripts"),
    function(x, transcripts, ...) 
        standardGeneric("pmapToTranscripts")
)

setGeneric("mapFromTranscripts", signature=c("x", "transcripts"),
    function(x, transcripts, ...) 
        standardGeneric("mapFromTranscripts")
)

setGeneric("pmapFromTranscripts", signature=c("x", "transcripts"),
    function(x, transcripts, ...) 
        standardGeneric("pmapFromTranscripts")
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
### mapToTranscripts() and mapFromTranscripts()
###

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

.mapFromTranscripts <- function(x, transcripts, hits, ignore.strand)
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

        ## non-hits
        seqname <- as.character(runValue(seqnames(transcripts)[txHits]))
        if (any(skip <- is.na(s) | is.na(e) | mismatch)) {
            s[skip] <- 0L 
            e[skip] <- -1L
            seqname[skip] <- "UNMAPPED"
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
    function(x, transcripts, ignore.strand=TRUE, ...)
    {
        grl <- relist(transcripts, PartitioningByEnd(seq_along(transcripts), 
                      names=names(transcripts)))
        mapToTranscripts(x, grl, ignore.strand, ...)
    }
)

setMethod("mapFromTranscripts", c("GenomicRanges", "GenomicRanges"), 
    function(x, transcripts, ignore.strand=TRUE, ...)
    {
        grl <- relist(transcripts, PartitioningByEnd(seq_along(transcripts), 
                      names=names(transcripts)))
        mapFromTranscripts(x, grl, ignore.strand, ...)
    }
)

setMethod("mapToTranscripts", c("GenomicRanges", "GRangesList"), 
    function(x, transcripts, ignore.strand=TRUE, ...) 
    {
        if (!length(x) && !length(transcripts))
            return(GRanges(xHits=integer(), transcriptsHits=integer()))

        if (!ignore.strand)
            if (!all(elementLengths(runLength(strand(transcripts))) == 1))
                stop(paste0("when ignore.strand=TRUE all inner list ",
                            "elements of 'transcripts' must be the ",
                            "same strand"))
        if (is.null(names(transcripts)))
            stop ("'transcripts' must have names")
        ## order within list elements by strand
        transcripts <- 
            .orderElementsByTranscription(transcripts, ignore.strand)
        ## findOverlaps determines pairs
        hits <- findOverlaps(x, unlist(transcripts, use.names=FALSE), 
                             type="within", ignore.strand=ignore.strand)
        .mapToTranscripts(x, transcripts, hits, ignore.strand)
    }
)

setMethod("mapFromTranscripts", c("GenomicRanges", "GRangesList"), 
    function(x, transcripts, ignore.strand=TRUE, ...) 
    {
        if (!length(x) && !length(transcripts))
            return(GRanges(xHits=integer(), transcriptsHits=integer()))

        if (!ignore.strand)
            if (!all(elementLengths(runLength(strand(transcripts))) == 1))
                stop(paste0("when ignore.strand=TRUE all inner list ",
                            "elements of 'transcripts' must be the ",
                            "same strand"))
        xNames <- names(x) 
        txNames <- names(transcripts)
        if (is.null(xNames) || is.null(txNames))
            stop ("both 'x' and 'transcripts' must have names")

        ## order within list elements by strand
        transcripts <- 
        .orderElementsByTranscription(transcripts, ignore.strand)
        ## name matching determines pairs
        match0 <- match(txNames, txNames)
        match1 <- match(xNames, txNames)
        group0 <- splitAsList(seq_along(txNames), match0)
        group1 <- group0[match(na.omit(match1), names(group0))]
        xHits <- rep(which(!is.na(match1)), elementLengths(group1))
        txHits <- unlist(group1, use.names=FALSE)
        if (!length(xHits <- na.omit(xHits)))
            stop ("none of 'names(x)' are in 'names(transcripts)'")

        hits <- Hits(xHits, txHits, length(x), length(transcripts))
        map <- .mapFromTranscripts(x, transcripts, hits, ignore.strand) 
        ## remove zero-width ranges 
        if (any(hasWidth <- width(map) != 0L))
            map <- map[hasWidth]
        map
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pmapToTranscripts() and pmapFromTranscripts()
###

setMethod("pmapToTranscripts", c("Ranges", "GenomicRanges"),
    function(x, transcripts, ...) 
    { 
        if (length(x) != length(transcripts))
            stop("'x' and 'transcripts' must have the same length")
        gr <- GRanges(seqnames(transcripts), x)
        ranges(pmapToTranscripts(gr, transcripts, TRUE))
    }
)

setMethod("pmapFromTranscripts", c("Ranges", "GenomicRanges"),
    function(x, transcripts, ...) 
    { 
        if (length(x) != length(transcripts))
            stop("'x' and 'transcripts' must have the same length")
        gr <- GRanges(seqnames(transcripts), x)
        ranges(pmapFromTranscripts(gr, transcripts, TRUE))
    }
)

setMethod("pmapToTranscripts", c("GenomicRanges", "GenomicRanges"), 
    function(x, transcripts, ignore.strand=TRUE, ...) 
    {
        grl <- splitAsList(transcripts, seq_along(transcripts))
        pmapToTranscripts(x, grl, ignore.strand, ...)
    }
)

setMethod("pmapFromTranscripts", c("GenomicRanges", "GenomicRanges"), 
    function(x, transcripts, ignore.strand=TRUE, ...) 
    {
        grl <- splitAsList(transcripts, seq_along(transcripts))
        pmapFromTranscripts(x, grl, ignore.strand, ...)
    }
)

setMethod("pmapToTranscripts", c("GenomicRanges", "GRangesList"), 
    function(x, transcripts, ignore.strand=TRUE, ...) 
    {
        if (!length(x) && !length(transcripts))
            return(GRanges(rep("UNMAPPED", length(x)), 
                           IRanges(rep(0, length(x)), width=0, names=names(x))))

        if (is.null(names(transcripts)))
            stop ("'transcripts' must have names")
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

        ## map i-th elements
        hits <- findOverlaps(x, unlist(transcripts, use.names=FALSE), 
                             type="within", ignore.strand=ignore.strand)
        ith <- queryHits(hits) == togroup(transcripts)[subjectHits(hits)]
        map <- .mapToTranscripts(x, transcripts, hits[ith], ignore.strand)

        ## non-hits
        if (length(x) != length(map)) {
            s <- rep(0L, length(x))
            e <- rep(-1L, length(x))
            seqname <- Rle("UNMAPPED", length(x))
            strands <- Rle("*", length(x))
            xHits <- mcols(map)$xHits
            e[xHits] <- end(map)
            s[xHits] <- start(map)
            seqname[xHits] <- seqnames(map)
            strands[xHits] <- strand(map)
            GRanges(seqname, IRanges(s, e, names=names(x)), strands)
        } else {
            mcols(map) <- NULL
            map
        }
    }
)

setMethod("pmapFromTranscripts", c("GenomicRanges", "GRangesList"), 
    function(x, transcripts, ignore.strand=TRUE, ...) 
    {
        if (!length(x) && !length(transcripts))
            return(GRanges(rep("UNMAPPED", length(x)), 
                           IRanges(rep(0, length(x)), width=0, names=names(x))))

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

        ## map i-th elements
        hits <- Hits(seq_along(x), seq_along(x), length(x), length(x))
        .mapFromTranscripts(x, transcripts, hits, ignore.strand)
    }
)
