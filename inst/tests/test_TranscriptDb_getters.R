###

test_getTranscripts <- function()
{
    txdb <- loadFeatures(system.file("extdata", "HG18test.sqlite", 
                                      package="GenomicFeatures"))
    suppressMessages(library(IRanges))
    ranges <- IRanges(start=c(1000), end=c(8000))
    chromosome <- "chr1"
    strand <-"+"
    wantRanges = IRanges(start = as.numeric(c(5658,4268,6720,6720)),
                         end   = as.numeric(c(7231,6628,7614,7924)))
    want<-RangedData(ranges     = wantRanges,
                     strand     = as.factor(c("-","-","-","-")),
                     space      = as.character(c("chr1","chr1","chr1","chr1")),
                     transcript = c("uc001aag.1","uc009vis.1",
                                    "uc009vjc.1","uc009vjd.1") )

    checkEquals(want, GenomicFeatures:::.getTranscripts(txdb,
                                                        ranges,
                                                        "chr1",
                                                        "-",
                                                        rangeRestr="both"))
}


test_getTranscripts <- function()
{
    txdb <- loadFeatures(system.file("extdata", "HG18test.sqlite", 
                                      package="GenomicFeatures"))
    suppressMessages(library(IRanges))
    rd <- RangedData(ranges = IRanges(start=c(1000), end=c(8000)),
                     strand = "-",
                     space  = "chr1")

    
    wantRanges = IRanges(start = as.numeric(c(5658,4268,6720,6720)),
                         end   = as.numeric(c(7231,6628,7614,7924)))
    want<-RangedData(ranges     = wantRanges,
                     strand     = as.factor(c("-","-","-","-")),
                     space      = as.character(c("chr1","chr1","chr1","chr1")),
                     transcript = c("uc001aag.1","uc009vis.1",
                                    "uc009vjc.1","uc009vjd.1") )

    checkEquals(want, getTranscripts(rd, txdb, rangeRestr="both") )
}

