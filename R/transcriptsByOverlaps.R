###

setGeneric("transcriptsByOverlaps", signature="x",
    function(x, ranges, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"), ...)
        standardGeneric("transcriptsByOverlaps")
)

setMethod("transcriptsByOverlaps", "TxDb",
    function(x, ranges, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             columns = c("tx_id", "tx_name"))
        subsetByOverlaps(transcripts(x, columns = columns), ranges,
                         maxgap = maxgap, minoverlap = minoverlap,
                         type = match.arg(type))
)

setGeneric("exonsByOverlaps", signature="x",
    function(x, ranges, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"), ...)
        standardGeneric("exonsByOverlaps")
)

setMethod("exonsByOverlaps", "TxDb",
    function(x, ranges, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             columns = "exon_id")
        subsetByOverlaps(exons(x, columns = columns), ranges,
                         maxgap = maxgap, minoverlap = minoverlap,
                         type = match.arg(type))
)

setGeneric("cdsByOverlaps", signature="x",
    function(x, ranges, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"), ...)
        standardGeneric("cdsByOverlaps")
)

setMethod("cdsByOverlaps", "TxDb",
    function(x, ranges, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             columns = "cds_id")
        subsetByOverlaps(cds(x, columns = columns), ranges,
                         maxgap = maxgap, minoverlap = minoverlap,
                         type = match.arg(type))
)

