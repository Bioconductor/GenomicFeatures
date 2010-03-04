### =========================================================================
### extractTranscriptsFromGenome()
### -------------------------------------------------------------------------


### 'x' must be a GRangesList object containing exons grouped by transcripts.
### Unlike the "as.data.frame" method for GRangesList, this function produces
### a data frame with 1 row per list element (a GRanges) and the following
### columns:
###   name       -- the outer name of the list element;
###   chrom      -- the single seqnames value of x[[i]] (raise an error
###                 if it has more than 1 value);
###   strand     -- the single strand value of x[[i]] (raise an error
###                 if it has more than 1 value);
###   exonStarts -- the starts of x[[i]] packed together in a single
###                 comma-separated string;
###   exonsEnds  -- the ends of x[[i]] packed together in a single
###                 comma-separated string.
.exonsByTxAsCompactDataFrame <- function(x, reorder.exons.on.minus.strand=TRUE)
{
    splitter <- rep.int(seq_len(length(x)), elementLengths(x))
    chrom <- split(as.character(seqnames(x@unlistData)), splitter)
    chrom <- sapply(chrom,
        function(y) {
            if (!all(y == y[1L]))
                stop("transcripts with exons from mixed chromosomes ",
                     "are not supported yet")
            y[1L]
        })
    strand <- split(as.character(strand(x@unlistData)), splitter)
    strand <- sapply(strand,
        function(y) {
            if (!all(y == y[1L]))
                stop("transcripts with exons from mixed strands ",
                     "are not supported yet")
            y[1L]
        })
    exonStarts <- split(start(x@unlistData), splitter)
    exonStarts <- sapply(seq_len(length(x)),
        function(i) {
            y <- exonStarts[[i]]
            if (reorder.exons.on.minus.strand && strand[i] == "-")
                y <- rev(y)
            paste(y, collapse=",")
        })
    exonEnds <- split(end(x@unlistData), splitter)
    exonEnds <- sapply(seq_len(length(x)),
        function(i) {
            y <- exonEnds[[i]]
            if (reorder.exons.on.minus.strand && strand[i] == "-")
                y <- rev(y)
            paste(y, collapse=",")
        })
    data.frame(
        name=names(x),
        chrom=unname(chrom),
        strand=unname(strand),
        exonStarts=unname(exonStarts),
        exonEnds=unname(exonEnds),
        stringsAsFactors=FALSE
    )   
}

### Typical use:
###   library(BSgenome.Hsapiens.UCSC.hg18)  # load the genome
###   library(GenomicFeatures.Hsapiens.UCSC.hg18)  # load the gene table
###   ## Takes about 30 sec.
###   transcripts <- extractTranscriptsFromGenome(Hsapiens, geneHuman())
### TODO: Improve implementation: ideally, we shouldn't need to use
###       intermediate representation as character vector.
extractTranscriptsFromGenome <- function(genome, txdb)
{
    if (!is(genome, "BSgenome"))
        stop("'genome' must be a BSgenome object")
    if (is(txdb, "TranscriptDb")) {
        exons_by_tx <- exonsBy(txdb, by="tx")
        genes <- .exonsByTxAsCompactDataFrame(exons_by_tx,
                     reorder.exons.on.minus.strand=FALSE)
        reorder.exons <- FALSE
    } else if (is.data.frame(txdb)) {
        REQUIRED_COLS <- c("name", "chrom", "strand", "exonStarts", "exonEnds")
        if (!all(REQUIRED_COLS %in% names(txdb)))
            stop("'txdb' data frame must have columns: ",
                 paste(REQUIRED_COLS, collapse=", "))
        genes <- txdb
        reorder.exons <- TRUE
    } else {
        stop("'txdb' must be a TranscriptDb object or a data frame")
    }
    ## The 3 lists below have identical names (the REFSEQnames)
    REFSEQnames2strand <- split(genes$strand, genes$chrom, drop=TRUE)
    REFSEQnames2exonStarts <- split(genes$exonStarts, genes$chrom, drop=TRUE)
    REFSEQnames2exonEnds <- split(genes$exonEnds, genes$chrom, drop=TRUE)
    REFSEQnames <- names(REFSEQnames2strand)  # REFSEQnames has no duplicates
    extractTranscriptSeqsFromREFSEQ <- function(REFSEQname)
    {
        subject <- genome[[REFSEQname]]
        masks(subject) <- NULL
        REFSEQ_exonStarts <- REFSEQnames2exonStarts[[REFSEQname]]
        REFSEQ_exonEnds <- REFSEQnames2exonEnds[[REFSEQname]]
        REFSEQ_strand <- REFSEQnames2strand[[REFSEQname]]
        transcripts <- extractTranscripts(subject,
                           REFSEQ_exonStarts, REFSEQ_exonEnds,
                           REFSEQ_strand,
                           reorder.exons.on.minus.strand=reorder.exons)
        as.character(transcripts)
    }
    REFSEQnames2seqs <- lapply(REFSEQnames, extractTranscriptSeqsFromREFSEQ)
    ans <- unsplit(REFSEQnames2seqs, genes$chrom, drop=TRUE)
    names(ans) <- genes$name
    DNAStringSet(ans)
}

