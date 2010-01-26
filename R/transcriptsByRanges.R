.txdbByRanges <-
function(txdb, ranges, restrict, FUN, prefix, ...)
{
    ## check that txdb is a TranscriptDb object
    if (!is(txdb, "TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")

    ## check that ranges is a RangedData object
    if (!is(ranges, "RangedData"))
        stop("'ranges' must be a RangedData object")

    useStrand <- ("strand" %in% colnames(ranges))

    chromName <- paste(prefix, "_chrom", sep="")
    do.call(c,
            lapply(names(ranges),
                   function(chrom) {
                       query <- ranges[chrom]
                       subject <-
                         FUN(txdb, structure(list(chrom), names=chromName), ...)
                       overlaps <-
                         findOverlaps(ranges(query)[[1L]],
                                      ranges(subject)[[1L]],
                                      type = restrict)
                       hits <- subjectHits(overlaps)
                       if (useStrand) {
                           hits <-
                             hits[strand(subject)[hits] ==
                                  strand(query)[queryHits(overlaps)]]
                       }
                       subject[sort(unique(hits)), , drop=FALSE]
                   }))
}


transcriptsByRanges <-
function(txdb, ranges, restrict = c("any", "start", "end", "within", "equal"),
         columns = c("tx_id", "tx_name"))
{
    .txdbByRanges(txdb=txdb, ranges=ranges, restrict=match.arg(restrict),
                  FUN=transcripts, prefix="tx", columns=columns)
}

exonsByRanges <-
function(txdb, ranges, restrict = c("any", "start", "end", "within", "equal"))
{
    .txdbByRanges(txdb=txdb, ranges=ranges, restrict=match.arg(restrict),
                  FUN=exons, prefix="exon")
}

cdsByRanges <-
function(txdb, ranges, restrict = c("any", "start", "end", "within", "equal"))
{
    .txdbByRanges(txdb=txdb, ranges=ranges, restrict=match.arg(restrict),
                  FUN=cds, prefix="cds")
}
