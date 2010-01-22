test_transcripts <- function()
{    
    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))
    suppressMessages(library(IRanges))
    vals <- list(tx_chrom = c("chr1", "chr5"), tx_strand = "-")
    
    wantRanges = IRanges(start = c(rep(4269,4),rep(257875,3)),
                         end   = c(rep(6628,4),rep(271297,3)))
    want<-RangedData(ranges  = wantRanges,
                     strand  = factor(rep("-",7),
                                     levels=c("-","+","*")),
                     space   = c(rep("chr1",4),rep("chr5",3)),
                     tx_id   = c(rep(110,4),rep(69,3)),
                     tx_name = c(rep("uc009vis.1",4),rep("uc003jam.1",3)))

    checkEquals(want, transcripts(txdb, vals))
}



test_exons <- function()
{
    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))
    suppressMessages(library(IRanges))
    vals <- list(tx_chrom = c("chr1", "chr5"), tx_strand = "-")
    
    wantRanges = IRanges(start = c(4269,4833,5659,6470,257875,269844,271208),
                         end   = c(4692,4901,5805,6628,259073,269974,271297))

    
    want<-RangedData(ranges  = wantRanges,
                     strand  = factor(rep("-",7),
                                     levels=c("-","+","*")),
                     space   = c(rep("chr1",4),rep("chr5",3)),
                     exon_id   = c(5,6,7,8,77,78,79))

    checkEquals(want, exons(txdb, vals))    
}
