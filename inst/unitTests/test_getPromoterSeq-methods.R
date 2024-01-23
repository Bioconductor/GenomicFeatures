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

test_GRangesListBSgenomeHumanGetPromoterSeq <- function() {
    txdb <- restoreSeqlevels(TxDb.Hsapiens.UCSC.hg19.knownGene)  ## safety net
    genes <- c(e2f3, grb2)
    tx_by_gene <- transcriptsBy(txdb, by="gene")[genes]
    checkIdentical(names(tx_by_gene), genes)
    transcript_count <- length(unlist(tx_by_gene, use.names=FALSE))

    promoter_seqs <- getPromoterSeq(tx_by_gene, Hsapiens,
                                    upstream=10, downstream=0)
    checkTrue(validObject(promoter_seqs))
    checkTrue(is(promoter_seqs, "DNAStringSetList"))
    checkEquals(length(promoter_seqs), 2)
    checkIdentical(names(promoter_seqs), genes)
    checkIdentical(width(unlist(promoter_seqs, use.names=FALSE)),
                   rep.int(10L, transcript_count))

    terminator_seqs <- getTerminatorSeq(tx_by_gene, Hsapiens,
                                        upstream=10, downstream=0)
    checkTrue(validObject(terminator_seqs))
    checkTrue(is(terminator_seqs, "DNAStringSetList"))
    checkEquals(length(terminator_seqs), 2)
    checkIdentical(names(terminator_seqs), genes)
    checkIdentical(width(unlist(terminator_seqs, use.names=FALSE)),
                   rep.int(10L, transcript_count))
}

test_GRangesListBSgenomeFlyGetPromoterSeq <- function() {
     # two neighboring genes near beginning of chr3R, on opposite strands
     #  gene_id  flybase_id  symbol
     #    40524 FBgn0037215 CG12582
     #    40526 FBgn0037217 CG14636
     # in 2012, UCSC reported 4 total transcripts for these two genes
     # in 2013, 6.  there should be as many promoter_seqs as there
     # are transcripts, and they should each be of width
     # upstream + downstream.  it is risky to check for specific
     # sequence in the promoter_seqs since the annotation and sequence
     # may change

    txdb <- restoreSeqlevels(TxDb.Dmelanogaster.UCSC.dm3.ensGene)  ## safety net
    genes <- c("FBgn0037215", "FBgn0037217")
    tx_by_gene <- transcriptsBy(txdb, by="gene")[genes]
    checkIdentical(names(tx_by_gene), genes)
    transcript_count <- length(unlist(tx_by_gene, use.names=FALSE))

    promoter_seqs <- getPromoterSeq(tx_by_gene, Dmelanogaster,
                                    upstream=10, downstream=10)
    checkTrue(validObject(promoter_seqs))
    checkTrue(is(promoter_seqs, "DNAStringSetList"))
    checkEquals(length(promoter_seqs), 2)
    checkIdentical(names(promoter_seqs), genes)
    checkIdentical(width(unlist(promoter_seqs, use.names=FALSE)),
                   rep.int(20L, transcript_count))

    terminator_seqs <- getPromoterSeq(tx_by_gene, Dmelanogaster,
                                      upstream=10, downstream=10)
    checkTrue(validObject(terminator_seqs))
    checkTrue(is(terminator_seqs, "DNAStringSetList"))
    checkEquals(length(terminator_seqs), 2)
    checkIdentical(names(terminator_seqs), genes)
    checkIdentical(width(unlist(terminator_seqs, use.names=FALSE)),
                   rep.int(20L, transcript_count))
}

test_GRangesListFastaFlyGetPromoterSeq <- function() {
      # two neighboring genes near beginning of chr3R, on opposite strands
      #  gene_id  flybase_id  symbol   chr
      #    43766 FBgn0025740  plexB    4
      #    43769 FBgn0085432  pan      4

    txdb <- restoreSeqlevels(TxDb.Dmelanogaster.UCSC.dm3.ensGene)  ## safety net
    genes <- c("FBgn0025740", "FBgn0085432")
    tx_by_gene <- transcriptsBy(txdb, by="gene")[genes]
    checkIdentical(names(tx_by_gene), genes)
    transcript_count <- length(unlist(tx_by_gene, use.names=FALSE))
    fa_file <- FaFile(dm3_chr4())

    promoter_seqs <- getPromoterSeq(tx_by_gene, fa_file,
                                    upstream=10, downstream=10)
    checkTrue(validObject(promoter_seqs))
    checkTrue(is(promoter_seqs, "DNAStringSetList"))
    checkEquals(length(promoter_seqs), 2)
    checkIdentical(names(promoter_seqs), genes)
    checkIdentical(width(unlist(promoter_seqs, use.names=FALSE)),
                   rep.int(20L, transcript_count))
       # we are unable to check for specific DNA sequence, since
       # the UCSC annotation of these genes changes over time.

    terminator_seqs <- getPromoterSeq(tx_by_gene, fa_file,
                                      upstream=10, downstream=10)
    checkTrue(validObject(terminator_seqs))
    checkTrue(is(terminator_seqs, "DNAStringSetList"))
    checkEquals(length(terminator_seqs), 2)
    checkIdentical(names(terminator_seqs), genes)
    checkIdentical(width(unlist(terminator_seqs, use.names=FALSE)),
                   rep.int(20L, transcript_count))
}

test_GRangesBSgenomeHumanGetPromoterSeq <- function() {
    txdb <- restoreSeqlevels(TxDb.Hsapiens.UCSC.hg19.knownGene)  ## safety net
    e2f3_tx <- transcriptsBy(txdb, by="gene")[[e2f3]]
    #names(e2f3_tx) <- mcols(e2f3_tx)$tx_name
    transcript_count <- length(e2f3_tx)
    checkEquals(dim(mcols(e2f3_tx)), c(transcript_count, 2))
    checkIdentical(colnames(mcols(e2f3_tx)), c("tx_id", "tx_name"))

    promoter_seqs <- getPromoterSeq(e2f3_tx, Hsapiens,
                                    upstream=10, downstream=0)
    checkTrue(validObject(promoter_seqs))
    checkTrue(is(promoter_seqs, "DNAStringSet"))
    checkEquals(length(promoter_seqs), transcript_count)
    checkTrue(is.null(names(promoter_seqs)))
    checkIdentical(width(promoter_seqs), rep.int(10L, transcript_count))
      # should be one more column in the metadata than in the metadata
    checkEquals(dim(mcols(promoter_seqs)), c(transcript_count, 3))
    checkEquals(colnames(mcols(promoter_seqs)), c("tx_id", "tx_name", "geneID"))
       # the input, a GRanges, had no names -- which are the source
       # of geneID when the GRangesList version of this methods is called.
       # so ensure that this lack of information was passed along into the
       # metadata of the returned promoter_seqs
    checkTrue(all(is.na(mcols(promoter_seqs)$geneID)))
}

