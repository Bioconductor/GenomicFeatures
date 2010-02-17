.txdbByRanges <-
function(txdb, ranges, restrict, FUN, prefix, ...)
{
    ## check that txdb is a TranscriptDb object
    if (!is(txdb, "TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")

    ## check that ranges is a GenomicFeature object
    if (!is(ranges, "GenomicFeature"))
        stop("'ranges' must be a GenomicFeature object")

    useStrand <- !all(is.na(strand(ranges)) | strand(ranges) == "*")

    chromName <- paste(prefix, "_chrom", sep="")
    do.call(c,
            lapply(unique(seqnames(ranges)),
                   function(chrom) {
                       query <- seqselect(ranges, seqnames(ranges) == chrom)
                       subject <-
                         FUN(txdb, structure(list(chrom), names=chromName), ...)
                       overlaps <-
                         findOverlaps(ranges(query), ranges(subject),
                                      type = restrict)
                       hits <- subjectHits(overlaps)
                       if (useStrand) {
                           hits <-
                             seqselect(hits,
                                       strand(subject)[hits] ==
                                         strand(query)[queryHits(overlaps)])
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
