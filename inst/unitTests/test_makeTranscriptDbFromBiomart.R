###

test_makeTranscriptDbFromBiomart <- function()
{
    ## want
    txdb0_file <- system.file(
                      "extdata",
                      "Biomart_Ensembl_sample.sqlite",
                      package="GenomicFeatures")
    txdb0 <- loadDb(txdb0_file)

    ## get
    transcript_ids <- c(
        "ENST00000268655",
        "ENST00000313243",
        "ENST00000341724",
        "ENST00000400839",
        "ENST00000435657",
        "ENST00000478783"
    )
    txdb1 <- makeTranscriptDbFromBiomart(transcript_ids=transcript_ids)

    ## compare
    ok <- GenomicFeatures:::compareTranscriptDbs(txdb1, txdb0)
    checkTrue(ok)
}

