### =========================================================================
### extractTranscriptsFromGenome()
### -------------------------------------------------------------------------


.exonsByTx <- function(txdb, use.tx_id=FALSE)
{
    if (!isTRUEorFALSE(use.tx_id))
        stop("'use.tx_id' must be TRUE or FALSE")
    exbytx <- exonsBy(txdb, by="tx")
    if (use.tx_id)
        return(exbytx)
    tx <- transcripts(txdb)
    id2name <- elementMetadata(tx)[ , "tx_name"]
    names(id2name) <- as.character(elementMetadata(tx)[ , "tx_id"])
    names(exbytx) <- id2name[names(exbytx)]
    if (any(is.na(names(exbytx))))
        warning("some transcript names are NAs, use 'use.tx_id=TRUE' ",
                "to use their internal ids instead")
    exbytx
}

### 'x' must be a GRangesList object containing exons grouped by transcripts.
### Unlike the "as.data.frame" method for GRangesList, this function produces
### a transcript table in the UCSC format i.e. a data frame with 1 row per
### element in 'x' (elements in 'x' are GRanges objects) and the following
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
### TODO: Improve implementation (very slow right-now).
.exonsByTxAsUCSCTxTable <- function(x, reorder.exons.on.minus.strand=TRUE)
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

.extractTranscriptsFromGenomeAndUCSCTxTable <- function(genome, ucsc_txtable,
                                                    reorder.exons.on.minus.strand)
{
    ## The 3 lists below have identical shapes and names (names are the
    ## REFSEQnames).
    strand_list <- split(ucsc_txtable$strand, ucsc_txtable$chrom, drop=TRUE)
    exonStarts_list <- split(ucsc_txtable$exonStarts, ucsc_txtable$chrom, drop=TRUE)
    exonEnds_list <- split(ucsc_txtable$exonEnds, ucsc_txtable$chrom, drop=TRUE)
    REFSEQnames <- names(strand_list)  # REFSEQnames has no duplicates
    extractTranscriptsFromREFSEQ <- function(REFSEQname)
    {
        subject <- genome[[REFSEQname]]
        masks(subject) <- NULL
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
                                      ucsc_txtable$chrom)
    names(ans) <- ucsc_txtable$name
    ans
}

### Typical use:
###   library(BSgenome.Hsapiens.UCSC.hg18)  # load the genome
###   library(GenomicFeatures.Hsapiens.UCSC.hg18)  # load the gene table
###   ## Takes about 30 sec.
###   transcripts <- extractTranscriptsFromGenome(Hsapiens, geneHuman())
extractTranscriptsFromGenome <- function(genome, txdb, use.tx_id=FALSE)
{
    if (!is(genome, "BSgenome"))
        stop("'genome' must be a BSgenome object")
    if (is(txdb, "TranscriptDb"))
        txdb <- .exonsByTx(txdb, use.tx_id=use.tx_id)
    if (is(txdb, "GRangesList")) {
        ucsc_txtable <- .exonsByTxAsUCSCTxTable(txdb,
                       reorder.exons.on.minus.strand=FALSE)
        reorder.exons <- FALSE
    } else if (is.data.frame(txdb)) {
        REQUIRED_COLS <- c("name", "chrom", "strand", "exonStarts", "exonEnds")
        if (!all(REQUIRED_COLS %in% names(txdb)))
            stop("'txdb' data frame must have columns: ",
                 paste(REQUIRED_COLS, collapse=", "))
        ucsc_txtable <- txdb
        reorder.exons <- TRUE
    } else {
        stop("'txdb' must be a TranscriptDb object or a data frame")
    }
    .extractTranscriptsFromGenomeAndUCSCTxTable(genome, ucsc_txtable, reorder.exons)
}

