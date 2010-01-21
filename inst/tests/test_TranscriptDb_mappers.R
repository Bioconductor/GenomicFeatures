test_mapTranscripts <- function()
{
    txdb <- loadFeatures(system.file("extdata", "HG18test.sqlite",
                                      package="GenomicFeatures"))
    suppressMessages(library(IRanges))
    ranges <- IRanges(start = c(1,20,300,1115),end =c(22,41,321,1115+21))
    chrom <- c("chr1", "chr1", "chr2", "chr1")
    strand <- c("+", "-", "+", "+")

    want <- RangedData(ranges     = IRanges(start = rep(1115,2),
                                        end   = rep(1136,2) ),
                   strand     = rep("+",2),
                   GF_txId    = c("uc009vip.1","uc001aaa.2"),
                   txName     = c("uc009vip.1","uc001aaa.2"),
                   space      = rep("chr1",2) )
    
    checkEquals(want, GenomicFeatures:::.mapTranscripts(txdb, ranges,
                                                        chrom, strand))

}


test_mapTranscripts <- function()
{
    txdb <- loadFeatures(system.file("extdata", "HG18test.sqlite",
                                      package="GenomicFeatures"))
    suppressMessages(library(IRanges))
    ranges <- IRanges(start = c(1,20,300,1115),end =c(22,41,321,1115+21))
    chrom <- c("chr1", "chr1", "chr2", "chr1")
    strand <- c("+", "-", "+", "+")

    want <- RangedData(ranges     = IRanges(start = rep(1115,2),
                                        end   = rep(1136,2) ),
                   strand     = rep("+",2),
                   GF_txId    = c("uc009vip.1","uc001aaa.2"),
                   txName     = c("uc009vip.1","uc001aaa.2"),
                   space      = rep("chr1",2) )

    rd <- RangedData(ranges, space = chrom, strand = strand)
    checkEquals(want, bindTranscripts(txdb, rd))
}





