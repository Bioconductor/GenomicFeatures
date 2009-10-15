### =========================================================================
### getTranscriptSeq()
### -------------------------------------------------------------------------

.splitChromCoords <- function(x) as.integer(strsplit(x, ",", fixed=TRUE)[[1]])

### Typical use:
###   library(GenomicFeatures.Hsapiens.UCSC.hg18)
###   library(BSgenome.Hsapiens.UCSC.hg18)
###   tseqs <- getTranscriptSeq(geneHuman(), Hsapiens)  # takes about 15 min.!
### TODO: Make it 100x faster!
getTranscriptSeq <- function(genes, genome)
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
        REFSEQ_strand <- REFSEQnames2strand[[REFSEQname]]
        REFSEQ_exonStarts <- REFSEQnames2exonStarts[[REFSEQname]]
        REFSEQ_exonEnds <- REFSEQnames2exonEnds[[REFSEQname]]
        if (!(REFSEQname %in% seqnames(genome))) {
            warning("transcripts from unknown chromosome \"", REFSEQname, "\" are reported as NAs")
            seqs <- character(length(REFSEQ_strand))
            seqs[] <- NA_character_
            return(seqs)
        }
        subject <- genome[[REFSEQname]]
        masks(subject) <- NULL
        seqs <- sapply(seq_len(length(REFSEQ_strand)),
                       function(i)
                       {
                           exons <- Views(subject,
                                          start=.splitChromCoords(REFSEQ_exonStarts[i]),
                                          end=.splitChromCoords(REFSEQ_exonEnds[i]))
                           transcript <- unlist(as(exons, "XStringSet"))
                           if (REFSEQ_strand[i] == "-")
                               transcript <- reverseComplement(transcript)
                           as.character(transcript)
                       })
        seqs
    }
    REFSEQnames2seqs <- lapply(REFSEQnames, extractTranscriptSeqsFromREFSEQ)
    ans <- unsplit(REFSEQnames2seqs, genes$chrom, drop=TRUE)
    names(ans) <- genes$name
    ans
}

