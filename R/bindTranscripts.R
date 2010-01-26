.txdbBind <-
function(txdb, ranges, restrict, bind, prefix, FUN, ...)
{
    ## check that txdb is a TranscriptDb object
    if (!is(txdb, "TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")

    ## check that ranges is a RangedData object
    if (!is(ranges, "RangedData"))
        stop("'ranges' must be a RangedData object")

    useStrand <- ("strand" %in% colnames(ranges))

    idName <- paste(prefix, "_id", sep="")
    nameName <- ifelse(prefix == "tx", paste(prefix, "_name", sep=""), "")
    chromName <- paste(prefix, "_chrom", sep="")
    rangesName <- paste(prefix, "_ranges", sep="")
    strandName <- paste(prefix, "_strand", sep="")
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
                       hitMatrix <- as.matrix(overlaps)
                       hitMatrix <-
                         hitMatrix[order(hitMatrix[,1L], hitMatrix[,2L]),,
                                   drop=FALSE]
                       if (useStrand) {
                           hitMatrix <-
                             hitMatrix[strand(subject)[subjectHits(overlaps)] ==
                                       strand(query)[queryHits(overlaps)], ,
                                       drop=FALSE]
                       }
                       alignedSubject <- subject[hitMatrix[,2L],,drop=FALSE]
                       queryGroup <-
                         factor(as.character(hitMatrix[,1L]),
                                levels=as.character(seq_len(nrow(query))))
                       if (idName %in% bind) {
                           query[[idName]] <-
                             unname(.newListBySplit("CompressedIntegerList",
                                                    alignedSubject[[idName]],
                                                    queryGroup))
                       }
                       if (nameName %in% bind) {
                           query[[nameName]] <-
                             unname(.newListBySplit("CompressedCharacterList",
                                                    alignedSubject[[nameName]],
                                                    queryGroup))
                       }
                       if (chromName %in% bind) {
                           query[[chromName]] <-
                             unname(.newListBySplit("CompressedCharacterList",
                                                    space(alignedSubject),
                                                    queryGroup))
                       }
                       if (rangesName %in% bind) {
                           query[[chromName]] <-
                             unname(.newListBySplit("CompressedIRangesList",
                                                    ranges(alignedSubject)[[1L]],
                                                    queryGroup))
                       }
                       if (strandName %in% bind) {
                           query[[strandName]] <-
                             unname(.newListBySplit("CompressedCharacterList",
                                                    as.character(strand(alignedSubject)),
                                                    queryGroup))
                       }
                       query
                   }))
}


bindTranscripts <-
function(txdb, ranges, restrict = c("any", "start", "end", "within", "equal"),
         columns = c("tx_id", "tx_name"))
{
    .txdbBind(txdb=txdb, ranges=ranges, restrict=match.arg(restrict),
              FUN=transcripts, prefix="tx", bind=columns,
              columns=intersect(columns, c("tx_id", "tx_name")))
}

bindExons <-
function(txdb, ranges, restrict = c("any", "start", "end", "within", "equal"),
         columns = "exon_id")
{
    .txdbBind(txdb=txdb, ranges=ranges, restrict=match.arg(restrict),
              FUN=exons, prefix="exon", bind=columns)
}

bindCDS <-
function(txdb, ranges, restrict = c("any", "start", "end", "within", "equal"),
         columns = "cds_id")
{
    .txdbBind(txdb=txdb, ranges=ranges, restrict=match.arg(restrict),
              FUN=cds, prefix="cds", bind=columns)
}
