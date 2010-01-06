###

test_getTranscripts <- function() {
    txdb <- loadFeatures(system.file("extdata", "HG18test.sqlite", 
                                      package="GenomicFeatures"))
    ranges = IRanges(start=c(1000), end=c(8000))

    want <- data.frame(
                       tx_id_ = as.integer(c(6,112,122,123)),
                       tx_id = c("uc001aag.1","uc009vis.1","uc009vjc.1","uc009vjd.1"),
                       chromosome = as.character(c("chr1","chr1","chr1","chr1")),
                       strand = as.character(c("-","-","-","-")),
                       tx_start = as.numeric(c(5658,4268,6720,6720)),
                       tx_end = as.numeric(c(7231,6628,7614,7924)),
                       stringsAsFactors = FALSE)

    names(want) = gsub("(.+)_$", "_\\1", names(want))
    
    checkEquals(want, getTranscripts(txdb, ranges, "chr1", "-",
                                     rangeRestr="both") )
}
