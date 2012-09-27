### =========================================================================
### Map internal ids to external names for a given feature type.
### -------------------------------------------------------------------------


id2name <- function(txdb, feature.type=c("tx", "exon", "cds"))
{
    if(!is(txdb,"TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")
    feature.type <- match.arg(feature.type)
    tablename <- switch(feature.type, tx="transcript", exon="exon", cds="cds")
    what_cols <- paste0(c("_", ""), feature.type, c("_id", "_name"))
    SQL <- paste("SELECT",
                 paste(what_cols, collapse=", "),
                 "FROM", tablename)
    data <- dbEasyQuery(AnnotationDbi:::dbConn(txdb), SQL)
    ans <- data[[what_cols[2L]]]
    names(ans) <- as.character(data[[what_cols[1L]]])
    ans
}

