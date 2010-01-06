###

test_makeTranscriptDbFromUCSCTxTable <- function()
{
    ucsc_txtable_file <- system.file(
                             "extdata",
                             "UCSCknownGene_sample.rda",
                             package="GenomicFeatures"
                         )
    load(ucsc_txtable_file)
    txdb1 <- GenomicFeatures:::.makeTranscriptDbFromUCSCTxTable(
                 UCSCknownGene_sample, "Human", "knownGene")
    txdb0_file <- system.file(
                      "extdata",
                      "UCSCknownGene_sample.sqlite",
                      package="GenomicFeatures")
    txdb0 <- loadFeatures(txdb0_file)
    ok <- GenomicFeatures:::compareTranscriptDbs(txdb1, txdb0)
    checkEquals(ok, TRUE)
}

