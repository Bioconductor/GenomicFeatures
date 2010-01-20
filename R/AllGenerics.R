### Should we rename the first arg. 'txdb' or just 'x'?
### Also, should we change the signature to dispatch on the first arg. only?

setGeneric("getTranscripts", signature="ranges",
           function(ranges=NULL, ann, rangeRestr="either", expand=FALSE)
           standardGeneric("getTranscripts"))

setGeneric("bindTranscripts", signature="ranges",
           function(txdb, ranges=NULL, restrict="any",
                    columns=c("tx_id", "tx_name"))
           standardGeneric("bindTranscripts"))






setGeneric("getExons", signature="ranges",
           function(ranges=NULL, ann, rangeRestr="either", expand=FALSE)
           standardGeneric("getExons"))

setGeneric("bindExons", signature="ranges",
           function(txdb, ranges=NULL, restrict="any",
                    columns=c("exon_id", "exon_name"))
           standardGeneric("exonsByRanges"))
