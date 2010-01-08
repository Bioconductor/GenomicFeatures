###

test_makeTranscriptDbFromUCSCTxTable <- function()
{
    ## want
    txdb0_file <- system.file(
                      "extdata",
                      "UCSCknownGene_sample.sqlite",
                      package="GenomicFeatures")
    txdb0 <- loadFeatures(txdb0_file)

    ## get
    ucsc_txtable_file <- system.file(
                             "extdata",
                             "UCSCknownGene_sample.rda",
                             package="GenomicFeatures"
                         )
    load(ucsc_txtable_file)
    txdb1 <- GenomicFeatures:::.makeTranscriptDbFromUCSCTxTable(
                 UCSCknownGene_sample, "Human", "knownGene")

    ## compare
    ok <- GenomicFeatures:::compareTranscriptDbs(txdb1, txdb0)
    checkEquals(ok, TRUE)
}

test_makeTranscriptDbFromBiomart <- function()
{
    ## want
    txdb0_file <- system.file(
                      "extdata",
                      "Biomart_Ensembl_sample.sqlite",
                      package="GenomicFeatures")
    txdb0 <- loadFeatures(txdb0_file)

    ## get
    some_ids <- c(
         "ENST00000400839",
         "ENST00000400840",
         "ENST00000478783",
         "ENST00000435657",
         "ENST00000268655",
         "ENST00000313243",
         "ENST00000341724"
    )
    txdb1 <- makeTranscriptDbFromBiomart(ensembl_transcript_ids=some_ids)

    ## compare
    ok <- GenomicFeatures:::compareTranscriptDbs(txdb1, txdb0)
    checkEquals(ok, TRUE)
}

