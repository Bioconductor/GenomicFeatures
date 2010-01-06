###

test_makeIdsForUniqueDataFrameRows <- function()
{
    x <- data.frame(
             chrom=c("chr2", "chr2", "chr2", "chr2", "chr1",
                     "chr2", "chr2", "chr1", "chr3", "chr1"),
             strand=c("+", "-", "-", "+", "+", "+", "+", "-", "-", "+"),
             start=c(5, 2, 2, 5, 4, 5, 5, 4, 2, 1),
             end=c(15, 12, 12, 15, 14, 13, 15, 14, 12, 11)
         )
    y <- unique(x)[GenomicFeatures:::makeIdsForUniqueDataFrameRows(x), ]
    row.names(y) <- NULL
    checkEquals(y, x)
}

