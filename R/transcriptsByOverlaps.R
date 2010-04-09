transcriptsByOverlaps <-
function(txdb, ranges, maxgap = 0L, minoverlap = 1L,
         type = c("any", "start", "end"), columns = c("tx_id", "tx_name"))
{
    subsetByOverlaps(transcripts(txdb, columns = columns), ranges,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type))
}

exonsByOverlaps <-
function(txdb, ranges, maxgap = 0L, minoverlap = 1L,
         type = c("any", "start", "end"))
{
    subsetByOverlaps(exons(txdb), ranges,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type))
}

cdsByOverlaps <-
function(txdb, ranges, maxgap = 0L, minoverlap = 1L,
         type = c("any", "start", "end"))
{
    subsetByOverlaps(cds(txdb), ranges,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type))
}
