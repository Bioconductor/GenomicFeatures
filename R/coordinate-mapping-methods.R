### =========================================================================
### mapToTranscriptome() and pmapToTranscriptome() methods
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generics 
###

setGeneric("mapToTranscriptome", signature=c("x", "transcriptome"),
    function(x, transcriptome, reverse=FALSE, ...) 
        standardGeneric("mapToTranscriptome")
)

setGeneric("pmapToTranscriptome", signature=c("x", "transcriptome"),
    function(x, transcriptome, reverse=FALSE, ...) 
        standardGeneric("pmapToTranscriptome")
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
### mapToTranscriptome()
###

### reverse = FALSE
.mapToTranscriptome <- function(x, transcriptome, hits, ignore.strand) 
{
    flat <- unlist(transcriptome)
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
        if (length(flat) > length(transcriptome)) {
            shifted <- .listCumsumShifted(width(transcriptome))
            xrange <- shift(xrange, shifted[txHits])
        }
        ## seqnames come from 'transcriptome'
        GRanges(rep(names(flat), elementLengths(flat))[txHits],
                xrange, strand(flat)[txHits],
                DataFrame(xHits, 
                          transcriptomeHits=togroup(transcriptome)[txHits])) 
    } else {
        ans <- GRanges()
        mcols(ans) <- DataFrame(xHits=integer(), transcriptomeHits=integer())
        ans
    }
}

### reverse = TRUE 
.rmapToTranscriptome <- function(x, transcriptome, hits, ignore.strand)
{
    if (length(hits)) {
        xHits <- queryHits(hits)
        txHits <- subjectHits(hits)
        ## strand mismatch
        if (!ignore.strand) {
            txstrand <- unlist(runValue(strand(transcriptome)), use.names=FALSE)
            mismatch <- (as.numeric(strand(x)[xHits]) + 
                         as.numeric(txstrand[txHits])) == 3L 
            strand <- as.character(txstrand)[txHits]
        } else {
            mismatch <- logical(length(xHits))
            strand <- rep("+", length(xHits))
        }
        xStart <- as.list(start(x)[xHits])
        xEnd <- as.list(end(x)[xHits])
        txStart <- as.list(start(transcriptome)[txHits])
        txEnd <- as.list(end(transcriptome)[txHits])
        s <- unlist(transcriptLocs2refLocs(xStart, txStart, txEnd, strand,
                                           FALSE, FALSE), use.names=FALSE)
        e <- unlist(transcriptLocs2refLocs(xEnd, txStart, txEnd, strand,
                                           FALSE, FALSE), use.names=FALSE)

        ## non-hits are zero width with seqname "unmmapped"
        seqname <- 
            as.character(runValue(seqnames(transcriptome)[txHits]))
        if (any(skip <- is.na(s) | is.na(e) | mismatch)) {
            s[skip] <- 1L 
            e[skip] <- 0L
            seqname[skip] <- "unmapped"
        }
        GRanges(Rle(seqname), IRanges(s, e, names=names(x)[xHits]), 
                strand=strand, DataFrame(xHits, transcriptomeHits=txHits))
    } else {
        ans <- GRanges()
        mcols(ans) <- DataFrame(xHits=integer(), transcriptomeHits=integer())
        ans
    }
}

setMethod("mapToTranscriptome", c("GenomicRanges", "GenomicRanges"), 
    function(x, transcriptome, reverse=FALSE, ignore.strand=TRUE, ...)
    {
        grl <- relist(transcriptome, PartitioningByEnd(seq_along(transcriptome), 
                      names=names(transcriptome)))
        if (reverse)
            mapToTranscriptome(x, grl, TRUE, ignore.strand, ...)
        else 
            mapToTranscriptome(x, grl, FALSE, ignore.strand, ...)
    }
)

setMethod("mapToTranscriptome", c("GenomicRanges", "GRangesList"), 
    function(x, transcriptome, reverse=FALSE, ignore.strand=TRUE, ...) 
    {
        if (length(x) && length(transcriptome)) {
            if (!ignore.strand)
                if (!all(elementLengths(runLength(strand(transcriptome))) == 1))
                    stop(paste0("when ignore.strand=TRUE all inner list ",
                                "elements of 'transcriptome' must be the ",
                                "same strand"))

            ## order within list elements by strand
            transcriptome <- 
                .orderElementsByTranscription(transcriptome, ignore.strand)
            if (reverse) {
                ## name matching determines pairs
                if (is.null(xNames <- names(x)) || 
                    is.null(transcriptomeNames <- names(transcriptome)))
                    stop ("both 'x' and 'transcriptome' must have names")
                match0 <- match(transcriptomeNames, transcriptomeNames)
                match1 <- match(xNames, transcriptomeNames)
                group0 <- splitAsList(seq_along(transcriptomeNames), match0)
                group1 <- group0[match(na.omit(match1), names(group0))]
                xHits <- rep(which(!is.na(match1)), elementLengths(group1))
                txHits <- unlist(group1, use.names=FALSE)
                if (!length(xHits <- na.omit(xHits)))
                    stop ("none of 'names(x)' are in 'names(transcriptome)'")

                hits <- Hits(xHits, txHits, length(x), length(transcriptome))
                map <- .rmapToTranscriptome(x, transcriptome, 
                                            hits, ignore.strand) 
                ## remove zero-width ranges 
                if (any(hasWidth <- width(map) != 0L))
                    map <- map[hasWidth]
                map
            } else {
                if (is.null(names(transcriptome)))
                    stop ("'transcriptome' must have names")
                ## findOverlaps determines pairs
                hits <- findOverlaps(x, unlist(transcriptome, use.names=FALSE), 
                                     type="within", ignore.strand=ignore.strand)
                .mapToTranscriptome(x, transcriptome, hits, ignore.strand)

            }
        } else {
            ans <- GRanges()
            mcols(ans) <- DataFrame(xHits=integer(), 
                                    transcriptomeHits=integer())
            ans
        }
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pmapToTranscriptome()
###

setMethod("pmapToTranscriptome", c("Ranges", "GenomicRanges"),
    function(x, transcriptome, reverse=FALSE, ...) 
    { 
        gr <- GRanges(seqnames(transcriptome), x)
        if (reverse)
            ranges(pmapToTranscriptome(gr, transcriptome, TRUE, TRUE))
        else 
            ranges(pmapToTranscriptome(gr, transcriptome, FALSE, TRUE))
    }
)

setMethod("pmapToTranscriptome", c("GenomicRanges", "GenomicRanges"), 
    function(x, transcriptome, reverse=FALSE,  ignore.strand=TRUE, ...) 
    {
        grl <- splitAsList(transcriptome, seq_along(transcriptome))
        if (reverse)
            pmapToTranscriptome(x, grl, TRUE, ignore.strand, ...)
        else 
            pmapToTranscriptome(x, grl, FALSE, ignore.strand, ...)
    }
)

setMethod("pmapToTranscriptome", c("GenomicRanges", "GRangesList"), 
    function(x, transcriptome, reverse=FALSE, ignore.strand=TRUE, ...) 
    {
        if (length(x) && length(transcriptome)) {
            if (length(x) != length(transcriptome))
                stop("'x' and 'transcriptome' must have the same length")
            if (!ignore.strand)
                if (!all(elementLengths(runLength(strand(transcriptome))) == 1))
                    stop(paste0("when ignore.strand=TRUE all inner list ",
                                "elements of 'transcriptome' must have the ",
                                "same strand"))
            ## order within list elements
            transcriptome <- 
                .orderElementsByTranscription(transcriptome, ignore.strand)

            if (reverse) {
                ## i-th element matching determines pairs
                hits <- Hits(seq_along(x), seq_along(x), length(x), length(x))
                .rmapToTranscriptome(x, transcriptome, hits, ignore.strand)
            } else {
                if (is.null(names(transcriptome)))
                    stop ("'transcriptome' must have names")
                ## i-th element matching determines pairs
                hits <- findOverlaps(x, unlist(transcriptome, use.names=FALSE), 
                                     type="within", ignore.strand=ignore.strand)
                ith <- 
                    queryHits(hits) == togroup(transcriptome)[subjectHits(hits)]
                map <- .mapToTranscriptome(x, transcriptome, 
                                           hits[ith], ignore.strand)

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

