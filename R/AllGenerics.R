### Should we rename the first arg. 'txdb' or just 'x'?
### Also, should we change the signature to dispatch on the first arg. only?

setGeneric("getTranscripts", signature="ranges",
           function(ranges=NULL, ann, rangeRestr="either", expand=FALSE)
           standardGeneric("getTranscripts"))

setGeneric("transcriptsByRanges", signature="ranges",
           function(txdb, ranges=NULL, restrict="any",
                    columns=c("tx_id", "tx_name"))
           standardGeneric("transcriptsByRanges"))






setGeneric("getExons", signature="ranges",
           function(ranges=NULL, ann, rangeRestr="either", expand=FALSE)
           standardGeneric("getExons"))

setGeneric("exonsByRanges", signature="ranges",
           function(txdb, ranges=NULL, restrict="any",
                    columns=c("exon_id", "exon_name"))
           standardGeneric("exonsByRanges"))
