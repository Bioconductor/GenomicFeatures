test__mapTranscripts <- function()
{
    txdb <- loadFeatures(system.file("extdata", "HG18test.sqlite",
                                      package="GenomicFeatures"))
    suppressMessages(library(IRanges))
    ranges <- IRanges(start = c(1,20,300,1115),end =c(22,41,321,1115+21))
    chromosome <- c("chr1", "chr1", "chr2", "chr1")
    strand <- c("+", "-", "+", "+")

    want <- RangedData(ranges     = IRanges(start = rep(1115,2),
                                        end   = rep(1136,2) ),
                   strand     = rep("+",2),
                   GF_txId    = as.integer(c(109,1)),
                   txId       = c("uc009vip.1","uc001aaa.2"),
                   space      = rep("chr1",2) )
    
    checkEquals(want, GenomicFeatures:::.mapTranscripts(txdb, ranges,
                                                        chromosome, strand))

}


test_mapTranscripts <- function()
{
    txdb <- loadFeatures(system.file("extdata", "HG18test.sqlite",
                                      package="GenomicFeatures"))
    suppressMessages(library(IRanges))
    ranges <- IRanges(start = c(1,20,300,1115),end =c(22,41,321,1115+21))
    chromosome <- c("chr1", "chr1", "chr2", "chr1")
    strand <- c("+", "-", "+", "+")

    want <- RangedData(ranges     = IRanges(start = rep(1115,2),
                                        end   = rep(1136,2) ),
                   strand     = rep("+",2),
                   GF_txId    = as.integer(c(109,1)),
                   txId       = c("uc009vip.1","uc001aaa.2"),
                   space      = rep("chr1",2) )

    rd <- RangedData(ranges, space = chromosome, strand = strand)
    checkEquals(want, mapTranscripts(rd, txdb))
}





