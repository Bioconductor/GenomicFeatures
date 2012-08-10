### =========================================================================
### Tools for extracting transcript sequences
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### transcriptLocs2refLocs()
###

transcriptLocs2refLocs <- function(tlocs, exonStarts=list(), exonEnds=list(),
                                   strand=character(0),
                                   decreasing.rank.on.minus.strand=FALSE)
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
                            decreasing.rank.on.minus.strand)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extractTranscripts()
###

extractTranscripts <- function(x, exonStarts=list(), exonEnds=list(),
                               strand=character(0),
                               decreasing.rank.on.minus.strand=FALSE)
{
    if (!is(x, "DNAString")) {
        if (!is(x, "MaskedDNAString"))
            stop("'x' must be a DNAString object")
        masks(x) <- NULL
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
    if (length(exonStarts) != length(strand)
     || length(exonEnds) != length(strand))
        stop("'exonStarts', 'exonEnds' and 'strand' must have the same length")
    if (!isTRUEorFALSE(decreasing.rank.on.minus.strand))
        stop("'decreasing.rank.on.minus.strand' must be TRUE or FALSE")
    lkup <- Biostrings:::getDNAComplementLookup()
    GenomicRanges:::unsafe.extractTranscripts("DNAStringSet", x,
                            exonStarts, exonEnds, strand,
                            decreasing.rank.on.minus.strand, lkup)
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
    exon_rank <- elementMetadata(x@unlistData)$exon_rank
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

.extractTranscriptsFromGenomeAndUCSCTxList <- function(genome, ucsc_txlist,
                                                decreasing.rank.on.minus.strand)
{
    ## The 3 lists below have identical shapes and names (names are the
    ## REFSEQnames).
    strand_list <- split(ucsc_txlist$strand, ucsc_txlist$chrom, drop=TRUE)
    exonStarts_list <- split(ucsc_txlist$exonStarts, ucsc_txlist$chrom, drop=TRUE)
    exonEnds_list <- split(ucsc_txlist$exonEnds, ucsc_txlist$chrom, drop=TRUE)
    REFSEQnames <- names(strand_list)  # REFSEQnames has no duplicates
    extractTranscriptsFromREFSEQ <- function(REFSEQname)
    {
        ## Load the subject.
        if (REFSEQname %in% seqnames(genome)) {
            subject <- genome[[REFSEQname]]
            masks(subject) <- NULL
        } else {
            regex <- paste0("^", REFSEQname, "$")
            subject <- getSeq(genome, regex, as.character=FALSE)
        }
        exonStarts <- exonStarts_list[[REFSEQname]]
        exonEnds <- exonEnds_list[[REFSEQname]]
        strand <- strand_list[[REFSEQname]]
        extractTranscripts(subject,
            exonStarts, exonEnds, strand,
            decreasing.rank.on.minus.strand=decreasing.rank.on.minus.strand)
    }
    ## Loop over the names of the reference sequences and extract the
    ## transcripts.
    dnaset_list <- lapply(REFSEQnames, extractTranscriptsFromREFSEQ)
    ans <- unsplit_list_of_XVectorList("DNAStringSet", dnaset_list,
                                       ucsc_txlist$chrom)
    names(ans) <- ucsc_txlist$name
    ans
}

extractTranscriptsFromGenome <- function(genome, txdb,
                                         decreasing.rank.on.minus.strand=FALSE,
                                         use.names=TRUE)
{
    if (!is(genome, "BSgenome"))
        stop("'genome' must be a BSgenome object")
    if (!isTRUEorFALSE(decreasing.rank.on.minus.strand))
        stop("'decreasing.rank.on.minus.strand' must be TRUE or FALSE")
    if (is.data.frame(txdb)) {
        ucsc_txlist <- .makeUCSCTxListFromUCSCTxTable(txdb)
    } else {
        if (is(txdb, "TranscriptDb")) {
            if (decreasing.rank.on.minus.strand)
                stop("'decreasing.rank.on.minus.strand' must be FALSE ",
                     "when 'txdb' is a TranscriptDb object")
            txdb <- exonsBy(txdb, by="tx", use.names=use.names)
        } else if (!is(txdb, "GRangesList"))
            stop("'txdb' must be a TranscriptDb object, a GRangesList ",
                 "object, or a data frame")
        ucsc_txlist <- .makeUCSCTxListFromGRangesList(txdb)
    }
    .extractTranscriptsFromGenomeAndUCSCTxList(genome, ucsc_txlist,
                                               decreasing.rank.on.minus.strand)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### sortExonsByRank()
###
### FIXME: Current implementation is a quick-and-dirty one that leverages the
### work done in .makeUCSCTxListFromGRangesList(). Therefore it looses the
### elementMetadata. We need something better, faster, and documented (and it
### should probably go somewhere else).
### TODO: Also, maybe exonsBy(... , by="tx") should get a
### 'by.decreasing.rank.on.minus.strand' arg or something like that.
###

sortExonsByRank <- function(x, decreasing.rank.on.minus.strand=FALSE)
{
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
                        strand=unlisted_strand)
    seqlevels(unlisted) <- seqlevels(x)
    seqinfo(unlisted) <- seqinfo(x)
    relist(unlisted, x)
}

