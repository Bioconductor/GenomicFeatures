test_transcriptLengths <- function()
{
    gff <- system.file("extdata", "GFF3_files",
                       "ITAG4.1_gene_models.subset.gff",
                       package="GenomicFeatures")
    txdb <- makeTxDbFromGFF(gff)
    txlens <- transcriptLengths(txdb, with.cds_len=TRUE,
                                      with.utr5_len=TRUE,
                                      with.utr3_len=TRUE)

    checkIdentical(class(txlens), "data.frame")
    checkIdentical(dim(txlens), c(10L, 8L))

    expected_colnames <- c("tx_id", "tx_name", "gene_id", "nexon", "tx_len",
                           "cds_len", "utr5_len", "utr3_len")
    checkIdentical(colnames(txlens), expected_colnames)

    checkIdentical(txlens$tx_len,
                   txlens$cds_len + txlens$utr5_len + txlens$utr3_len)

    expected_nexon <- c(2L, 1L, 1L, 3L, 1L, 2L, 3L, 4L, 2L, 2L)
    checkIdentical(txlens$nexon, expected_nexon)

    expected_tx_len <- c(721L, 813L, 319L, 1046L, 264L,
                         819L, 802L, 453L, 249L, 889L)
    checkIdentical(txlens$tx_len, expected_tx_len)

    expected_cds_len <- c(471L, 663L, 276L, 918L, 264L,
                          261L, 369L, 453L, 249L, 516L)
    checkIdentical(txlens$cds_len, expected_cds_len)

    expected_utr5_len <- c(250L, 0L, 0L, 89L, 0L, 558L, 0L, 0L, 0L, 373L)
    checkIdentical(txlens$utr5_len, expected_utr5_len)

    expected_utr3_len <- c(0L, 150L, 43L, 39L, 0L, 0L, 433L, 0L, 0L, 0L)
    checkIdentical(txlens$utr3_len, expected_utr3_len)
}

