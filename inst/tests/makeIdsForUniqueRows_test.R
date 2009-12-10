makeIdsForUniqueRows <- GenomicFeatures:::makeIdsForUniqueRows

test_makeIdForUniqueRows <- function() {
    df <- data.frame(chrom=c("chr2", "chr2", "chr2", "chr2", "chr1", "chr2",
                     "chr2", "chr1", "chr3", "chr1"),
                     strand=c("+", "-", "-", "+", "+", "+", "+", "-", "-",
                     "+"),
                     start=c(5, 2, 2, 5, 4, 5, 5, 4, 2, 1),
                     end=c(15, 12, 12, 15, 14, 13, 15, 14, 12, 11)) 

    want <- c(6,4,4,6,3,5,6,1,7,2)
  
    checkEquals(want, makeIdsForUniqueRows(df))
}
