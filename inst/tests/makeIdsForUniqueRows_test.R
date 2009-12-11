set.seed(0xab56)

makeIdsForUniqueRows <- GenomicFeatures:::makeIdsForUniqueRows


test_makeIdForUniqueRows_extraCols <- function()
{
    df <- data.frame(foo=c("a", "b", "c"), 
                     chrom=c("chr2", "chr1", "chr2"),
                     strand=c("+", "-", "+"),
                     sna=LETTERS[4:6],
                     exonStart=c(5, 2, 5),
                     exonEnd=c(15, 12, 15),
                     bar=1:3)
    want <- c(2, 1, 2)
    checkEquals(want, makeIdsForUniqueRows(df))

    ord <- c(3, 2, 1)
    checkEquals(want[ord], makeIdsForUniqueRows(df[ord, ]))

    ord <- c(2, 3, 1)
    checkEquals(want[ord], makeIdsForUniqueRows(df[ord, ]))

    ord <- c(2, 1, 3)
    checkEquals(want[ord], makeIdsForUniqueRows(df[ord, ]))
}





test_makeIdForUniqueRows_simple <- function()
{
    df <- data.frame(chrom=c("chr2", "chr1", "chr2"),
                     strand=c("+", "-", "+"),
                     exonStart=c(5, 2, 5),
                     exonEnd=c(15, 12, 15))
    want <- c(2, 1, 2)
    checkEquals(want, makeIdsForUniqueRows(df))

    ord <- c(3, 2, 1)
    checkEquals(want[ord], makeIdsForUniqueRows(df[ord, ]))

    ord <- c(2, 3, 1)
    checkEquals(want[ord], makeIdsForUniqueRows(df[ord, ]))

    ord <- c(2, 1, 3)
    checkEquals(want[ord], makeIdsForUniqueRows(df[ord, ]))
}

test_makeIdForUniqueRows <- function() {
    df <- data.frame(chrom=c("chr2", "chr2", "chr2", "chr2", "chr1", "chr2",
                     "chr2", "chr1", "chr3", "chr1"),
                     strand=c("+", "-", "-", "+", "+", "+", "+", "-", "-",
                     "+"),
                     exonStart=c(5, 2, 2, 5, 4, 5, 5, 4, 2, 1),
                     exonEnd=c(15, 12, 12, 15, 14, 13, 15, 14, 12, 11))

    ord <- do.call(order, df)
    df <- df[ord, ]
    isDup <- duplicated(df)
    want <- vector(mode="numeric", length=nrow(df))
    id <- 0
    for (i in seq_len(nrow(df))) {
        if (!isDup[i]) id <- id + 1
        want[i] = id
    }

    checkEquals(want, makeIdsForUniqueRows(df))

    for (i in 1:50) {
        ord <- sample(1:nrow(df))
        checkEquals(want[ord], makeIdsForUniqueRows(df[ord, ]))
    }
}
