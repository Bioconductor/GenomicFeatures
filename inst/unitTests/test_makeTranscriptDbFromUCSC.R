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
                 "hg18", "knownGene", "Entrez Gene ID", FALSE,
                 DEFAULT_CIRC_SEQS)

    ## compare
    ok <- GenomicFeatures:::compareTranscriptDbs(txdb1, txdb0)
    checkTrue(ok)
}

