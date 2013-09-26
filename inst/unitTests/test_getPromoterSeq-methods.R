library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(Rsamtools)
library(pasillaBamSubset)

e2f3 <- "1871"   # human gene on the plus strand, chr6
grb2 <- "2885"   # human gene on the minus strand, chr17

# a note on method: when the promoter sequence is 20 bases or more in length,
# uscs blat will find these sequences, and a quick visual inspection of the
# accompanying genome browser view at the right level of zoom, will
# confirm that the per-transcript sequences is indeed correct.
# there are a few tests of shorter sequences below as well, which
# I checked in the genome browser, but this required a little more effort
# than the length 20, blat approach.

testGRangesListBSgenomeHumanGetPromoterSeq <- function() {
    genes <- c(e2f3, grb2)
    TxDb.Hsapiens.UCSC.hg19.knownGene <- restoreSeqlevels(TxDb.Hsapiens.UCSC.hg19.knownGene)  ## safety net
    transcriptCoordsByGene.GRangesList <-
      transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="gene") [genes]
    transcript.count <- length(unlist(transcriptCoordsByGene.GRangesList))

    checkEquals(names(transcriptCoordsByGene.GRangesList), genes)
    promoter.seqs <- getPromoterSeq(transcriptCoordsByGene.GRangesList,
                                    Hsapiens, upstream=10, downstream=0)
    checkTrue(is(promoter.seqs, "DNAStringSetList"))
    checkEquals(length(promoter.seqs), 2)
    checkEquals(names(promoter.seqs), genes)
    checkEquals(width(unlist(promoter.seqs)), rep(10, transcript.count))
}

testGRangesListBSgenomeFlyGetPromoterSeq <- function() {
     # two neighboring genes near beginning of chr3R, on opposite strands
     #  gene_id  flybase_id  symbol
     #    40524 FBgn0037215 CG12582
     #    40526 FBgn0037217 CG14636
     # in 2012, UCSC reported 4 total transcripts for these two genes
     # in 2013, 6.  there should be as many promoter.seqs as there
     # are transcripts, and they should each be of width
     # upstream + downstream.  it is risky to check for specific
     # sequence in the promoter.seqs since the annotation and sequence
     # may change

    genes <- c("FBgn0037215", "FBgn0037217")
    transcriptCoordsByGene.GRangesList <-
       transcriptsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, by="gene") [genes]
  
    transcript.count <- length(unlist(transcriptCoordsByGene.GRangesList))

    promoter.seqs <- getPromoterSeq(transcriptCoordsByGene.GRangesList,
                                    Dmelanogaster, upstream=10, downstream=10)
    checkTrue(is(promoter.seqs, "DNAStringSetList"))
    checkEquals(length(promoter.seqs), 2)
    checkEquals(names(promoter.seqs), genes)
  
    checkEquals(width(unlist(promoter.seqs)), rep(20, transcript.count))
}

testGRangesListFastaFlyGetPromoterSeq <- function() {
      # two neighboring genes near beginning of chr3R, on opposite strands
      #  gene_id  flybase_id  symbol   chr
      #    43766 FBgn0025740  plexB    4
      #    43769 FBgn0085432  pan      4
    genes <- c("FBgn0025740", "FBgn0085432")
    transcriptCoordsByGene.GRangesList <-
       transcriptsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, by="gene") [genes]

    transcript.count <- length(unlist(transcriptCoordsByGene.GRangesList))

    fasta.file <- dm3_chr4 ()
    sequence.from.fasta <- open(FaFile(fasta.file))
    promoter.seqs <- getPromoterSeq(transcriptCoordsByGene.GRangesList,
                                    sequence.from.fasta, upstream=10,
                                    downstream=10)
    checkTrue(is(promoter.seqs, "DNAStringSetList"))
    checkEquals(length(promoter.seqs), 2)
    checkEquals(names(promoter.seqs), genes)
  
    checkEquals(width(unlist(promoter.seqs)), rep(20, transcript.count))
       # we are unable to check for specific DNA sequence, since
       # the UCSC annotation of these genes changes over time.
}

testGRangesBSgenomeHumanGetPromoterSeq <- function() {
    TxDb.Hsapiens.UCSC.hg19.knownGene <- restoreSeqlevels(TxDb.Hsapiens.UCSC.hg19.knownGene)  ## safety net
    transcriptCoordsByGene.GRanges <-
      transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="gene") [[e2f3]]
    checkTrue(is(transcriptCoordsByGene.GRanges, "GRanges"))
       # would have names only if its a list:
    transcript.count <- length(transcriptCoordsByGene.GRanges)
    
    checkTrue(is.null(names(transcriptCoordsByGene.GRanges)))
    checkEquals(dim(mcols(transcriptCoordsByGene.GRanges)),
                c(transcript.count, 2))
    checkEquals(colnames(mcols(transcriptCoordsByGene.GRanges)),
                c("tx_id", "tx_name"))
    promoter.seqs <-
      getPromoterSeq(transcriptCoordsByGene.GRanges, Hsapiens,
                     upstream=10, downstream=0)
    checkTrue(is(promoter.seqs, "DNAStringSet"))
    checkEquals(length(promoter.seqs), transcript.count)
    checkTrue(is.null(names(promoter.seqs)))
    checkEquals(width(promoter.seqs), rep(10, transcript.count))
      # should be one more column in the metadata than in the metadata 
    checkEquals(dim(mcols(promoter.seqs)), c(transcript.count, 3))
    checkEquals(colnames(mcols(promoter.seqs)), c("tx_id", "tx_name", "geneID"))
       # the input, a GRanges, had no names -- which are the source
       # of geneID when the GRangesList version of this methods is called.
       # so ensure that this lack of information was passed along into the
       # metadata of the returned promoter.seqs
    checkTrue(all(is.na(mcols(promoter.seqs)$geneID)))
}

