### =========================================================================
### extractTranscriptsFromGenome()
### -------------------------------------------------------------------------

### Typical use:
###   library(BSgenome.Hsapiens.UCSC.hg18)  # load the genome
###   library(GenomicFeatures.Hsapiens.UCSC.hg18)  # load the gene table
###   ## Takes about 30 sec.
###   transcripts <- extractTranscriptsFromGenome(Hsapiens, geneHuman())
### TODO: Improve implementation: ideally, we shouldn't need to use
###       intermediate representation as character vector.
extractTranscriptsFromGenome <- function(genome, genes)
{
    if (!is(genome, "BSgenome"))
        stop("'genome' must be a BSgenome object")
    if (is.data.frame(genes)) {
        REQUIRED_COLS <- c("name", "chrom", "strand", "exonStarts", "exonEnds")
        if (!all(REQUIRED_COLS %in% names(genes)))
            stop("'genes' data frame must have columns: ",
                 paste(REQUIRED_COLS, collapse=", "))
        reorder.exons <- TRUE
    } else if (is(genes, "TranscriptDb")) {
        exons_by_tx <- exonsBy(genes, by="tx")
        chrom <- sapply(exons_by_tx,
            function(tx) {
                ans <- seqnames(tx)
                if (length(runValue(ans)) != 1L)
                    stop("transcripts with exons from mixed chroms ",
                         "are not supported yet")
                as.character(runValue(ans))
            })
        strand <- sapply(exons_by_tx,
            function(tx) {
                ans <- strand(tx)
                if (length(runValue(ans)) != 1L)
                    stop("transcripts with exons from mixed strands ",
                         "are not supported yet")
                as.character(runValue(ans))
            })
        exonStarts <- sapply(exons_by_tx,
            function(tx) paste(start(tx), collapse=","))
        exonEnds <- sapply(exons_by_tx,
            function(tx) paste(end(tx), collapse=","))
        genes <- data.frame(
            name=names(exons_by_tx),
            chrom=unname(chrom),
            strand=unname(strand),
            exonStarts=unname(exonStarts),
            exonEnds=unname(exonEnds),
            stringsAsFactors=FALSE
        )   
        reorder.exons <- FALSE
    } else {
        stop("'genes' must be a data frame or a TranscriptDb object")
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

