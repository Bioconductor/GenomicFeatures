###

test_getTranscripts <- function()
{
    txdb <- loadFeatures(system.file("extdata", "HG18test.sqlite", 
                                      package="GenomicFeatures"))
    suppressMessages(library(IRanges))
    ranges <- IRanges(start=1000, end=8000)
    chrom <- "chr1"
    strand <-"+"
    wantRanges = IRanges(start = c(5658,4268,6720,6720),
                         end   = c(7231,6628,7614,7924))
    wantIds <- c("uc001aag.1","uc009vis.1", "uc009vjc.1","uc009vjd.1")
    want<-RangedData(ranges = wantRanges,
                     strand = factor(rep("-",4),
                                     levels=c("-","+","*")),
                     space  = rep("chr1",4),
                     GF_txId = wantIds,
                     txName = wantIds)

    checkEquals(want, GenomicFeatures:::.mapTranscripts(txdb,
                                                        ranges,
                                                        "chr1",
                                                        "-",
                                                        type="any"))
}


test_getTranscripts <- function()
{
    txdb <- loadFeatures(system.file("extdata", "HG18test.sqlite", 
                                      package="GenomicFeatures"))
    suppressMessages(library(IRanges))
    rd <- RangedData(ranges = IRanges(start=1000, end=8000),
                     space  = "chr1")    
    wantRanges = IRanges(start = c(1115,5658,1115,4268,6720,6720),
                         end   = c(4121,7231,4272,6628,7614,7924))
    wantIds <- c("uc001aaa.2","uc001aag.1","uc009vip.1",
                 "uc009vis.1","uc009vjc.1","uc009vjd.1")
    want<-RangedData(ranges = wantRanges,
                     strand = factor(c("+","-","+","-","-","-"),
                                     levels=c("-","+","*")),
                     space  = rep("chr1",6),
                     GF_txId = wantIds,
                     txName = wantIds)

    checkEquals(want, transcriptsByRanges(txdb, rd, restrict="any") )
}

