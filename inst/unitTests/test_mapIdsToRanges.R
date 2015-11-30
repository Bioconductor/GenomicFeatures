txdb <- makeTxDbFromGRanges(
    readRDS(system.file(package = "GenomicFeatures", "extdata", "sample_ranges.rds")))

test_mapIdsToRanges_improper_inputs <- function()
{
    checkException(mapIdsToRanges(txdb, vals = ""),
        "must be a named list")

    checkException(mapIdsToRanges(txdb, vals = "ENST000000271582"),
        "must be a named list")

    checkException(mapIdsToRanges(txdb, vals = list("ENST000000271582")),
        "must be a named list")

    checkException(mapIdsToRanges(txdb, vals = list(tx_name = "ENST000000271582"), column = 1),
        "'columns' must be 'NULL' or a character vector")
}

test_mapIdsToRanges_same_order <- function()
{
    vals <- list(tx_name = c("ENST00000371582", "ENST00000371588", "ENST00000494752", "ENST00000614008", "ENST00000496771"))
    res <- mapIdsToRanges(txdb, vals = vals, type = "tx")
    checkEquals(names(res), vals[[1]])

    # shuffle the order and make sure it remains equivalent
    for (i in seq_len(10)) {
        vals$tx_name <- sample(vals$tx_name)
        res <- mapIdsToRanges(txdb, vals = vals, type = "tx")
        checkEquals(names(res), vals[[1]])
    }
}

test_mapIdsToRanges_missing_results <- function()
{
    vals <- list(tx_name = c("ENST00000371582", "NOT_FOUND", "ENST00000494752"))
    res <- mapIdsToRanges(txdb, vals = vals, type = "tx")
    checkEquals(names(res), vals$tx_name)

    # shuffle the order and make sure it remains equivalent
    for (i in seq_len(10)) {
        vals$tx_name <- sample(vals$tx_name)
        res <- mapIdsToRanges(txdb, vals = vals, type = "tx")
        checkEquals(names(res), vals$tx_name)
    }
}

test_mapIdsToRanges_duplicate_ranges <- function()
{
    # both of these transcripts are from the same gene
    vals <- list(tx_name = c("ENST00000371582", "ENST00000494752"))
    res <- mapIdsToRanges(txdb, vals = vals, type = "gene")

    #names match input
    checkEquals(names(res), vals[[1]])

    # but values are the same
    checkTrue(all.equal(res[[1]], res[[2]], check.attributes = FALSE))
}

test_mapRangesToIds_empty <- function()
{
    checkException(mapRangesToIds(txdb, NULL), list())

    checkException(mapRangesToIds(txdb, list()), list())
}

test_mapRangesToIds_matches <- function()
{
    vals <- list(tx_name = c("ENST00000371582", "ENST00000371588", "ENST00000494752", "ENST00000614008", "ENST00000496771"))
    res <- mapIdsToRanges(txdb, vals = vals, type = "gene")

    res2 <- mapRangesToIds(txdb, res, "tx")
    checkTrue(vals$tx_name[1] %in% res2[[1]]$tx_name)

    checkTrue(vals$tx_name[2] %in% res2[[2]]$tx_name)
}
