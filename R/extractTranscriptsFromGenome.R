### =========================================================================
### extractTranscriptsFromGenome()
### -------------------------------------------------------------------------

### Typical use:
###   library(GenomicFeatures.Hsapiens.UCSC.hg18)
###   library(BSgenome.Hsapiens.UCSC.hg18)
###   ## Takes < 30 sec.
###   tseqs <- extractTranscriptsFromGenome(geneHuman(), Hsapiens)
### TODO: Improve implementation: ideally, we shouldn't need to use
###       intermediate representation as character vector.
extractTranscriptsFromGenome <- function(genes, genome)
{
    if (!is.data.frame(genes))
        stop("'genes' must be a data.frame")
    if (!is(genome, "BSgenome"))
        stop("'genome' must be a BSgenome object")
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
        exonStarts <- strsplit(REFSEQ_exonStarts, ",", fixed=TRUE)
        exonStarts <- lapply(exonStarts, as.integer)
        exonEnds <- strsplit(REFSEQ_exonEnds, ",", fixed=TRUE)
        exonEnds <- lapply(exonEnds, as.integer)
        transcripts <- extractTranscripts(subject, exonStarts, exonEnds,
                           REFSEQ_strand, reorder.exons.on.minus.strand=TRUE)
        as.character(transcripts)
    }
    REFSEQnames2seqs <- lapply(REFSEQnames, extractTranscriptSeqsFromREFSEQ)
    ans <- unsplit(REFSEQnames2seqs, genes$chrom, drop=TRUE)
    names(ans) <- genes$name
    DNAStringSet(ans)
}

