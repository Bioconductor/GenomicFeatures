###

setMethod("transcriptsByOverlaps", "TranscriptDb",
    function(x, ranges, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             columns = c("tx_id", "tx_name"))
        subsetByOverlaps(transcripts(x, columns = columns), ranges,
                         maxgap = maxgap, minoverlap = minoverlap,
                         type = match.arg(type))
)

setMethod("exonsByOverlaps", "TranscriptDb",
    function(x, ranges, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             columns = "exon_id")
        subsetByOverlaps(exons(x, columns = columns), ranges,
                         maxgap = maxgap, minoverlap = minoverlap,
                         type = match.arg(type))
)

setMethod("cdsByOverlaps", "TranscriptDb",
    function(x, ranges, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             columns = "cds_id")
        subsetByOverlaps(cds(x, columns = columns), ranges,
                         maxgap = maxgap, minoverlap = minoverlap,
                         type = match.arg(type))
)
