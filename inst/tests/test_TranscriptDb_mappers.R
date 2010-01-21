test_mapTranscripts <- function()
{
    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                      package="GenomicFeatures"))
    suppressMessages(library(IRanges))
    ranges <- IRanges(start = c(1,20,300,31000),end =c(22,41,321,37000))
    chrom <- c("chr1", "chr1", "chr2", "chr2")
    strand <- c("-", "-", "+", "-")

    want <- RangedData(ranges     = IRanges(start = rep(31000,2),
                                            end   = rep(37000,2) ),
                   strand   = factor(rep("-",2),
                                     levels=c("-","+","*")),
                   tx_id    = rep(37L,2),
                   tx_name  = rep("uc002qvt.1",2),
                   space    = rep("chr2",2) )
    
    checkEquals(want, GenomicFeatures:::.mapTranscripts(txdb, ranges,
                                                        chrom, strand))

}


test_mapTranscripts <- function()
{
    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                      package="GenomicFeatures"))
    suppressMessages(library(IRanges))
    ranges <- IRanges(start = c(1,20,300,31000),end =c(22,41,321,37000))
    chrom <- c("chr1", "chr1", "chr2", "chr2")
    strand <- c("-", "-", "+", "-")

    want <- RangedData(ranges     = IRanges(start = rep(31000,2),
                                            end   = rep(37000,2) ),
                   strand   = factor(rep("-",2),
                                     levels=c("-","+","*")),
                   tx_id    = rep(37L,2),
                   tx_name  = rep("uc002qvt.1",2),
                   space    = rep("chr2",2) )
    
    rd <- RangedData(ranges, space = chrom, strand = strand)
    checkEquals(want, bindTranscripts(txdb, rd))
}


