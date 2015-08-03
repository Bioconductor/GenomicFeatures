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
### This function returns a GRangesList with sorted elements. It 
### differs from sort() in that "-" strand elements are returned 
### highest value to lowest.
.orderElementsByTranscription <- function(x) {
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

    neg <- strand(gr)[pstart] == "-"
    ord <- S4Vectors:::mseq(ifelse(neg, pend, pstart),
                            ifelse(neg, pstart, pend))
    res <- relist(gr[ord], x)
    res@unlistData$unordered <- original[idx[ord]] 
    res
}

### 'x' is an IntegerList or NumericList
### This function returns a numeric vector of cumulative sums within list 
### elements.
.listCumsumShifted <- function(x) {
    cs <- unlist(cumsum(x), use.names=FALSE)
    shifted <- c(0L, head(cs, -1))
    shifted[start(PartitioningByWidth(elementLengths(x)))] <- 0L
    shifted
}

.recycling <- function(xlen, txlen)
{
    if (xlen != txlen) {
        if (txlen == 1L)
            TRUE
        else
            stop(paste0("recycling is supported when length(transcripts) == 1; ",
                        "otherwise lengths of 'x' and 'transcripts' ",
                        "must match"))
    } else FALSE
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToTranscripts()
###

.mapToTranscripts <- function(x, transcripts, hits, ignore.strand) 
{
    flat <- unlist(transcripts, use.names=FALSE)
    if (length(hits)) {
        xHits <- queryHits(hits)
        txHits <- subjectHits(hits)
        xrange <- ranges(x)[xHits]
        bounds <- ranges(flat)[txHits]
        ## location wrt to start of individual list elements
        neg <- as.vector(strand(flat)[txHits] == "-")
        negstart <- end(bounds)[neg] - end(xrange)[neg] + 1L
        xrange[neg] <- IRanges(negstart, width=width(xrange)[neg])
        xrange[!neg] <- shift(xrange[!neg], - start(bounds)[!neg] + 1L)
        ## location wrt start of concatenated list elements
        if (length(flat) > length(transcripts)) {
            shifted <- .listCumsumShifted(width(transcripts))
            xrange <- shift(xrange, shifted[txHits])
        }
        ## seqnames come from 'transcripts'
        txGroup <- togroup(transcripts)[txHits]
        GRanges(names(transcripts)[txGroup], 
                xrange, strand(flat)[txHits],
                DataFrame(xHits, transcriptsHits=txGroup))
    } else {
        ans <- GRanges()
        mcols(ans) <- DataFrame(xHits=integer(), transcriptsHits=integer())
        ans
    }
}

setMethod("mapToTranscripts", c("GenomicRanges", "GenomicRanges"), 
    function(x, transcripts, ignore.strand=FALSE, ...)
    {
        grl <- relist(transcripts, PartitioningByEnd(seq_along(transcripts), 
                      names=names(transcripts)))
        mapToTranscripts(x, grl, ignore.strand, ...)
    }
)

setMethod("mapToTranscripts", c("GenomicRanges", "GRangesList"), 
    function(x, transcripts, ignore.strand=FALSE, ...) 
    {
        if (!length(x) && !length(transcripts))
            return(GRanges(xHits=integer(), transcriptsHits=integer()))

        if (!ignore.strand)
            if (any(elementLengths(runValue(strand(transcripts))) > 1))
                stop(paste0("when ignore.strand=TRUE all inner list ",
                            "elements of 'transcripts' must be the ",
                            "same strand"))
        if (is.null(names(transcripts)))
            stop ("'transcripts' must have names")
        ## order within list elements by strand
        transcripts <- 
            .orderElementsByTranscription(transcripts)
        ## findOverlaps determines pairs
        hits <- findOverlaps(x, unlist(transcripts, use.names=FALSE), 
                             type="within", ignore.strand=ignore.strand)
        .mapToTranscripts(x, transcripts, hits, ignore.strand)
    }
)

setMethod("mapToTranscripts", c("ANY", "TxDb"), 
    function(x, transcripts, ignore.strand=FALSE, 
             extractor.fun = GenomicFeatures::transcripts, ...) 
    {
        if (!is.function(extractor.fun))
            stop("'extractor.fun' must be a function")
        group1 <- c("transcripts", "exons", "cds", "genes", "promoters",
                    "disjointExons", "microRNAs", "tRNAs") 
        group2 <- c("transcriptsBy", "exonsBy", "cdsBy", "intronsByTranscript", 
                   "fiveUTRsByTranscript", "threeUTRsByTranscript")

        fname <- extractor.fun@generic 
        if (fname %in% group1) {
            transcripts <- extractor.fun(transcripts, ...)
            if (is.null(names(transcripts)))
                names(transcripts) <- mcols(transcripts)[,1]
        } else if (fname %in% group2) {
            transcripts <- extractor.fun(transcripts, ...)
        } else {
            stop("invalid 'extractor.fun'")
        }
        mapToTranscripts(x, transcripts, ignore.strand, ...)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapFromTranscripts()
###

## use seqnames of 'x' and names of 'transcripts' for mapping

.mapFromTranscripts <- function(x, transcripts, hits, ignore.strand)
{
    if (length(hits)) {
        xHits <- queryHits(hits)
        txHits <- subjectHits(hits)
        ## strand mismatch
        txstrand <- unlist(runValue(strand(transcripts)), use.names=FALSE)
        strand <- tmpstrand <- as.character(txstrand)[txHits]
        if (ignore.strand || all(txstrand == "*")) {
            mismatch <- logical(length(xHits))
            tmpstrand <- rep("+", length(xHits)) ## for transcriptLocs2refLocs
        } else if (!ignore.strand) {
            mismatch <- (as.numeric(strand(x)[xHits]) + 
                         as.numeric(txstrand[txHits])) == 3L 
        }
        xStart <- as.list(start(x)[xHits])
        xEnd <- as.list(end(x)[xHits])
        txStart <- as.list(start(transcripts)[txHits])
        txEnd <- as.list(end(transcripts)[txHits])
        s <- unlist(transcriptLocs2refLocs(xStart, txStart, txEnd, tmpstrand,
                                           FALSE, FALSE), use.names=FALSE)
        e <- unlist(transcriptLocs2refLocs(xEnd, txStart, txEnd, tmpstrand,
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

setMethod("mapFromTranscripts", c("GenomicRanges", "GenomicRanges"), 
    function(x, transcripts, ignore.strand=FALSE, ...)
    {
        grl <- relist(transcripts, PartitioningByEnd(seq_along(transcripts), 
                      names=names(transcripts)))
        mapFromTranscripts(x, grl, ignore.strand, ...)
    }
)

setMethod("mapFromTranscripts", c("GenomicRanges", "GRangesList"), 
    function(x, transcripts, ignore.strand=FALSE, ...) 
    {
        if (!length(x) && !length(transcripts))
            return(GRanges(xHits=integer(), transcriptsHits=integer()))

        if (!ignore.strand)
            if (!all(elementLengths(runLength(strand(transcripts))) == 1))
                stop(paste0("when ignore.strand=TRUE all inner list ",
                            "elements of 'transcripts' must be the ",
                            "same strand"))
        xNames <- as.character(seqnames(x))
        txNames <- names(transcripts)
        if (is.null(txNames))
            stop ("'transcripts' must have names")

        ## order within list elements by strand
        transcripts <- 
        .orderElementsByTranscription(transcripts)
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
### pmapToTranscripts()
###

.pmapToTranscripts_ranges <- function(x, transcripts, ...)
{
    if (!length(x))
        return(IRanges())
    if (.recycling(length(x), length(transcripts)))
        transcripts <- rep(transcripts, length(x))
    gr <- GRanges(seqnames(transcripts), x)
    ranges(pmapToTranscripts(gr, transcripts, TRUE))
}

setMethod("pmapToTranscripts", c("Ranges", "GenomicRanges"),
    function(x, transcripts, ...) .pmapToTranscripts_ranges
)

setMethod("pmapToTranscripts", c("Ranges", "GRangesList"), 
    function(x, transcripts, ...) .pmapToTranscripts_ranges
)

setMethod("pmapToTranscripts", c("GenomicRanges", "GenomicRanges"), 
    function(x, transcripts, ignore.strand=FALSE, ...) 
    {
        grl <- splitAsList(transcripts, seq_along(transcripts))
        pmapToTranscripts(x, grl, ignore.strand, ...)
    }
)

setMethod("pmapToTranscripts", c("GenomicRanges", "GRangesList"), 
    function(x, transcripts, ignore.strand=FALSE, ...) 
    {
        if (!length(x))
            return(GRanges())
        if (.recycling(length(x), length(transcripts)))
            transcripts <- rep(transcripts, length(x))
        if (is.null(names(transcripts)))
            stop ("'transcripts' must have names")
        if (!ignore.strand)
            if (!all(elementLengths(runLength(strand(transcripts))) == 1))
                stop(paste0("when ignore.strand=TRUE all inner list ",
                            "elements of 'transcripts' must have the ",
                            "same strand"))
        ## order within list elements
        transcripts <- 
            .orderElementsByTranscription(transcripts)

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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pmapFromTranscripts()
###

.pmapFromTranscripts_ranges <- function(x, transcripts, ...) 
{
    if (!length(x))
        return(GRanges())
    if (.recycling(length(x), length(transcripts)))
        transcripts <- rep(transcripts, length(x))

    ## mapping
    xStart <- as.list(start(x))
    xEnd <- as.list(end(x))
    txStart <- as.list(start(transcripts))
    txEnd <- as.list(end(transcripts))
    tmpstrand <- rep("+", length(x)) ## for transcriptLocs2refLocs
    s <- unlist(transcriptLocs2refLocs(xStart, txStart, txEnd, tmpstrand,
                                       FALSE, FALSE), use.names=FALSE)
    e <- unlist(transcriptLocs2refLocs(xEnd, txStart, txEnd, tmpstrand,
                                       FALSE, FALSE), use.names=FALSE)

    ## non-hits
    if (any(skip <- is.na(s) | is.na(e))) {
        s[skip] <- 0L 
        e[skip] <- -1L
        seqname[skip] <- "UNMAPPED"
    }
    IRanges(s, e, names=names(x)) 
}

setMethod("pmapFromTranscripts", c("Ranges", "GenomicRanges"),
    function(x, transcripts, ...) 
        .pmapFromTranscripts_ranges
)

setMethod("pmapFromTranscripts", c("Ranges", "GRangesList"), 
    function(x, transcripts, ...) 
        .pmapFromTranscripts_ranges 
)

setMethod("pmapFromTranscripts", c("GenomicRanges", "GenomicRanges"), 
    function(x, transcripts, ignore.strand=FALSE, ...) 
    {
        grl <- splitAsList(transcripts, seq_along(transcripts))
        pmapFromTranscripts(x, grl, ignore.strand, ...)
    }
)

setMethod("pmapFromTranscripts", c("GenomicRanges", "GRangesList"), 
    function(x, transcripts, ignore.strand=FALSE, ...) 
    {
        if (!length(x))
            return(GRanges())
        if (.recycling(length(x), length(transcripts)))
            transcripts <- rep(transcripts, length(x))
        if (!ignore.strand)
            if (!all(elementLengths(runLength(strand(transcripts))) == 1))
                stop(paste0("when ignore.strand=TRUE all inner list ",
                            "elements of 'transcripts' must have the ",
                            "same strand"))
        ## order within list elements
        transcripts <- 
            .orderElementsByTranscription(transcripts)

        ## map i-th elements
        hits <- Hits(seq_along(x), seq_along(x), length(x), length(x))
        .mapFromTranscripts(x, transcripts, hits, ignore.strand)
    }
)
