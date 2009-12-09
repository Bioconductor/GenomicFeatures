convertExonsCommaSepFrame <- GenomicFeatures:::convertExonsCommaSepFrame

test_simpleConvert <- function() {
    eStarts <- c(paste(1:3, collapse=","),
                 100L,
                 paste(10:11, collapse=","))
    eEnds <- c(paste(1:3 + 10L, collapse=","),
               110L,
               paste(10:11 + 10L, collapse=","))

    df <- data.frame(apples = letters[1:3],
                     exonStarts = eStarts,
                     exonEnds = eEnds,
                     stringsAsFactors = FALSE)

    want <- data.frame(
                       apples = c(rep("a", 3L), "b", rep("c", 2L)),
                       exonStarts = c(1L, 2L, 3L, 100L, 10L, 11L),
                       exonEnds = c(1L, 2L, 3L, 100L, 10L, 11L) + 10L,
                       exonRank = as.integer(c(1, 2, 3, 1, 1, 2)),
                       stringsAsFactors = FALSE)

    checkEquals(want, convertExonsCommaSepFrame(df))
}
