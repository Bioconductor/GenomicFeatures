### =========================================================================
### extractTranscriptsFromGenome()
### -------------------------------------------------------------------------


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
.makeUCSCTxListFromGRangesList <- function(x, reorder.exons.on.minus.strand=TRUE)
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
                "inferring rank from exon position within GRanges")
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
    if (!is.null(exon_rank) || reorder.exons.on.minus.strand) {
        exonStarts <- sapply(seq_len(length(x)),
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
                if (reorder.exons.on.minus.strand && strand[i] == "-")
                    y <- rev(y)
                y
            })
        exonEnds <- sapply(seq_len(length(x)),
            function(i) {
                y <- exonEnds[[i]]
                if (!is.null(exon_rank)) {
                    perm <- exon_rank[[i]]
                    shift <- min(perm) - 1L
                    if (shift != 0L)
                        perm <- perm - shift
                    y[perm] <- y
                }
                if (reorder.exons.on.minus.strand && strand[i] == "-")
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
                                                  reorder.exons.on.minus.strand)
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
            regex <- paste("^", REFSEQname, "$", sep="")
            subject <- getSeq(genome, regex, as.character=FALSE)
        }
        exonStarts <- exonStarts_list[[REFSEQname]]
        exonEnds <- exonEnds_list[[REFSEQname]]
        strand <- strand_list[[REFSEQname]]
        extractTranscripts(subject,
            exonStarts, exonEnds, strand,
            reorder.exons.on.minus.strand=reorder.exons.on.minus.strand)
    }
    ## Loop over the names of the reference sequences and extract the
    ## transcripts.
    dnaset_list <- lapply(REFSEQnames, extractTranscriptsFromREFSEQ)
    ans <- unsplit.list.of.XStringSet("DNAStringSet", dnaset_list,
                                      ucsc_txlist$chrom)
    names(ans) <- ucsc_txlist$name
    ans
}

extractTranscriptsFromGenome <- function(genome, txdb, use.names=TRUE)
{
    if (!is(genome, "BSgenome"))
        stop("'genome' must be a BSgenome object")
    if (is.data.frame(txdb)) {
        ucsc_txlist <- .makeUCSCTxListFromUCSCTxTable(txdb)
        reorder.exons <- TRUE
    } else {
        if (is(txdb, "TranscriptDb")) {
            txdb <- exonsBy(txdb, by="tx", use.names=use.names)
        } else if (!is(txdb, "GRangesList"))
            stop("'txdb' must be a TranscriptDb object, a GRangesList ",
                 "object, or a data frame")
        ucsc_txlist <- .makeUCSCTxListFromGRangesList(txdb,
                           reorder.exons.on.minus.strand=FALSE)
        reorder.exons <- FALSE
    }
    .extractTranscriptsFromGenomeAndUCSCTxList(genome, ucsc_txlist, reorder.exons)
}

