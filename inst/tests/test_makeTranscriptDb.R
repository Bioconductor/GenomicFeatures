###

test_makeTranscriptDbFromUCSCTxTable <- function()
{
    ## want
    txdb0_file <- system.file(
                      "extdata",
                      "UCSC_knownGene_sample.sqlite",
                      package="GenomicFeatures"
                  )
    txdb0 <- loadFeatures(txdb0_file)

    ## get
    knownGene_sample_file <- system.file(
                                 "extdata",
                                 "UCSC_knownGene_sample.rda",
                                 package="GenomicFeatures"
                             )
    load(knownGene_sample_file)
    knownToLocusLink_sample_file <- system.file(
                                        "extdata",
                                        "UCSC_knownToLocusLink_sample.rda",
                                        package="GenomicFeatures"
                                    )
    load(knownToLocusLink_sample_file)
    genes <- data.frame(tx_name=UCSC_knownToLocusLink_sample$name,
                        gene_id=UCSC_knownToLocusLink_sample$value)
    txdb1 <- GenomicFeatures:::.makeTranscriptDbFromUCSCTxTable(
                 UCSC_knownGene_sample, genes,
                 "hg18", "knownGene", "Entrez Gene ID", FALSE)

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
    transcript_ids <- c(
         "ENST00000400839",
         "ENST00000400840",
         "ENST00000478783",
         "ENST00000435657",
         "ENST00000268655",
         "ENST00000313243",
         "ENST00000341724"
    )
    txdb1 <- makeTranscriptDbFromBiomart(transcript_ids=transcript_ids)

    ## compare
    ok <- GenomicFeatures:::compareTranscriptDbs(txdb1, txdb0)
    checkEquals(ok, TRUE)
}

test_makeTranscriptDb <- function()
{
    ## want
    txdb0_file <- system.file(
                      "extdata",
                      "Biomart_Ensembl_sample.sqlite",
                      package="GenomicFeatures")
    txdb0 <- loadFeatures(txdb0_file)

    ## get
    txdb1 <- do.call(makeTranscriptDb, as.list(txdb0))

    ## compare
    ok <- GenomicFeatures:::compareTranscriptDbs(txdb1, txdb0)
    checkEquals(ok, TRUE)
}

