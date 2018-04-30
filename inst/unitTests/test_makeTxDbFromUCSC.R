###

if (FALSE) {
test_makeTxDbFromUCSC <- function()
{
    txdb_file <- system.file("extdata",
                             "sacCer2_sgdGene_txdb.sqlite",
                             package="GenomicFeatures")
    target_txdb <- loadDb(txdb_file)

    current_txdb <- makeTxDbFromUCSC("sacCer2", "sgdGene")

    ok <- GenomicFeatures:::compareTxDbs(target_txdb, current_txdb)
    checkTrue(ok)
}
}
