### =========================================================================
### Map internal ids to external names for a given feature type.
### -------------------------------------------------------------------------


id2name <- function(txdb, feature.type=c("tx", "exon", "cds"))
{
    if (!is(txdb, "TxDb"))
        stop("'txdb' must be a TxDb object")
    feature.type <- match.arg(feature.type)
    table <- switch(feature.type, tx="transcript", exon="exon", cds="cds")
    columns <- TXDB_table_columns(table)[c("id", "name")]
    df <- TxDb_SELECT_from_INNER_JOIN(txdb, table, columns)
    ans <- df[[columns[2L]]]
    names(ans) <- as.character(df[[columns[1L]]])
    ans
}

