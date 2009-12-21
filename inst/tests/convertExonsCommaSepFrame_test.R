convertExonsCommaSepFrame <- GenomicFeatures:::convertExonsCommaSepFrame

test_simpleConvert <- function() {
    eStarts <- c(paste(1:3, collapse=","),
                 100L,
                 paste(10:11, collapse=","))
    eEnds <- c(paste(1:3 + 10L, collapse=","),
               110L,
               paste(10:11 + 10L, collapse=","))

    df <- data.frame(apples = letters[1:3],
                     exonStart = eStarts,
                     exonEnd = eEnds,
                     stringsAsFactors = FALSE)

    want <- data.frame(
                       apples = c(rep("a", 3L), "b", rep("c", 2L)),
                       exonStart = c(1L, 2L, 3L, 100L, 10L, 11L),
                       exonEnd = c(1L, 2L, 3L, 100L, 10L, 11L) + 10L,
                       exonRank = as.integer(c(1, 2, 3, 1, 1, 2)),
                       stringsAsFactors = FALSE)

    checkEquals(want, convertExonsCommaSepFrame(df))
}
