### =========================================================================
### Tools for extracting transcript sequences
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractTranscripts()
###

extractTranscripts <- function(x, exonStarts=list(), exonEnds=list(),
                               strand=character(0),
                               decreasing.rank.on.minus.strand=FALSE)
{
    .Defunct("extractTranscriptSeqs")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractTranscriptsFromGenome()
###

### Turns data frame 'ucsc_txtable' into a list where the "exonStarts" and
### "exonsEnds" columns (character vectors where each element is a
### comma-separated list of integers) are expanded into lists of integer
### vectors with no NAs.
.makeUCSCTxListFromUCSCTxTable <- function(ucsc_txtable)
{
    REQUIRED_COLS <- c("name", "chrom", "strand", "exonStarts", "exonEnds")
    if (!all(REQUIRED_COLS %in% names(ucsc_txtable)))
        stop("data frame must have columns: ",
             paste(REQUIRED_COLS, collapse=", "))
    list(name=ucsc_txtable$name,
         chrom=ucsc_txtable$chrom,
         strand=ucsc_txtable$strand,
         exonStarts=strsplitAsListOfIntegerVectors(ucsc_txtable$exonStarts),
         exonEnds=strsplitAsListOfIntegerVectors(ucsc_txtable$exonEnds))
}

### 'x' must be a GRangesList object containing exons (or CDSs) grouped by
### transcripts.
### Returns a named list with 5 elements of the same length (like the
### .makeUCSCTxListFromUCSCTxTable() function above). This common length is
### the length of 'x'. The elements are:
###   name       -- Character vector. 'names(x)'.
###   chrom      -- Character vector. 'chrom[i]' is the unique seqnames value
###                 found in 'x[[i]]' (an error is raised if 'x[[i]]' has more
###                 then 1 distinct seqnames value).
###   strand     -- Character vector with only possible values being "+" or
###                 "-". 'strand[i]' is the unique strand value found
###                 in 'x[[i]]' (an error is raised if 'x[[i]]' has more than
###                 1 distinct strand value).
###   exonStarts -- List of integer vectors with no NAs. Vector
###                 'exonStarts[[i]]' is made of the starts of 'x[[i]]',
###                 eventually ordered according to the exon rank information
###                 found in 'x[[i]]'.
###   exonsEnds  -- List of integer vectors with no NAs. Vector
###                 'exonEnds[[i]]' is made of the ends of 'x[[i]]',
###                 eventually ordered according to the exon rank information
###                 found in 'x[[i]]'.
.makeUCSCTxListFromGRangesList <- function(x,
                                      decreasing.rank.on.minus.strand=FALSE)
{
    f <- rep.int(seq_len(length(x)), elementLengths(x))
    ## Note that 'x@unlistData' is 50000x faster than
    ## 'unlist(x, use.names=FALSE)' and 3 million times faster
    ## than 'unname(unlist(x))'.
    chrom <- unname(split(as.character(seqnames(x@unlistData)), f))
    chrom <- sapply(chrom,
        function(y) {
            if (!all(y == y[1L]))
                stop("transcripts with exons from mixed chromosomes ",
                     "are not supported yet")
            y[1L]
        })
    strand <- unname(split(as.character(strand(x@unlistData)), f))
    strand <- sapply(strand,
        function(y) {
            if (!all(y == y[1L]))
                stop("transcripts with exons from mixed strands ",
                     "are not supported yet")
            y[1L]
        })
    exon_rank <- mcols(x@unlistData)$exon_rank
    if (is.null(exon_rank)) {
        warning("GRangesList object has no \"exon_rank\" column --> ",
                "inferring rank from exon position within each GRanges")
    } else {
        if (!is.numeric(exon_rank))
            stop("\"exon_rank\" column in GRangesList object is not numeric")
        if (!is.integer(exon_rank)) {
            warning("\"exon_rank\" column in GRangesList object is not integer")
            exon_rank <- as.integer(exon_rank)
        }
        exon_rank <- unname(split(exon_rank, f))
    }
    exonStarts <- unname(split(start(x@unlistData), f))
    exonEnds <- unname(split(end(x@unlistData), f))
    if (!is.null(exon_rank) || decreasing.rank.on.minus.strand) {
        exonStarts <- lapply(seq_len(length(x)),
            function(i) {
                y <- exonStarts[[i]]
                if (!is.null(exon_rank)) {
                    perm <- exon_rank[[i]]
                    ## When 'x' contains CDSs grouped by transcripts, some of
                    ## the lowest or/and highest exon ranks can be missing.
                    ## In that case, 'shift' will be != 0.
                    shift <- min(perm) - 1L
                    if (shift != 0L)
                        perm <- perm - shift
                    y[perm] <- y
                }
                if (decreasing.rank.on.minus.strand && strand[i] == "-")
                    y <- rev(y)
                y
            })
        exonEnds <- lapply(seq_len(length(x)),
            function(i) {
                y <- exonEnds[[i]]
                if (!is.null(exon_rank)) {
                    perm <- exon_rank[[i]]
                    shift <- min(perm) - 1L
                    if (shift != 0L)
                        perm <- perm - shift
                    y[perm] <- y
                }
                if (decreasing.rank.on.minus.strand && strand[i] == "-")
                    y <- rev(y)
                y
            })
    }
    list(name=names(x),
         chrom=chrom,
         strand=strand,
         exonStarts=exonStarts,
         exonEnds=exonEnds)
}

extractTranscriptsFromGenome <- function(genome, txdb,
                                         decreasing.rank.on.minus.strand=FALSE,
                                         use.names=TRUE)
{
    .Defunct("extractTranscriptSeqs")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sortExonsByRank()
###
### FIXME: Current implementation is a quick-and-dirty one that leverages the
### work done in .makeUCSCTxListFromGRangesList(). Therefore it looses the
### metadata columns. We need something better, faster, and documented (and it
### should probably go somewhere else).
### TODO: Also, maybe exonsBy(... , by="tx") should get a
### 'by.decreasing.rank.on.minus.strand' arg or something like that.
###

sortExonsByRank <- function(x, decreasing.rank.on.minus.strand=FALSE)
{
    .Deprecated()
    if (!is(x, "GRangesList"))
        stop("'x' must be a GRangesList object, ",
             "typically coming from exonsBy(... , by=\"tx\")")
    if (!isTRUEorFALSE(decreasing.rank.on.minus.strand))
        stop("'decreasing.rank.on.minus.strand' must be TRUE or FALSE")
    ucsc_txlist <- .makeUCSCTxListFromGRangesList(x,
        decreasing.rank.on.minus.strand=decreasing.rank.on.minus.strand)
    nexon <- elementLengths(ucsc_txlist$exonStarts)
    unlisted_seqnames <- Rle(ucsc_txlist$chrom, nexon)
    unlisted_strand <- Rle(strand(ucsc_txlist$strand), nexon)
    unlisted_start <- unlist(ucsc_txlist$exonStarts, use.names=FALSE)
    unlisted_end <- unlist(ucsc_txlist$exonEnds, use.names=FALSE)
    unlisted <- GRanges(seqnames=unlisted_seqnames,
                        ranges=IRanges(unlisted_start, unlisted_end),
                        strand=unlisted_strand,
                        seqinfo=seqinfo(x))
    relist(unlisted, x)
}
