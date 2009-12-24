setGeneric("getTranscripts", signature=c("ann", "ranges"),
           function(ann, ranges=NULL, chromosome=NULL, strand=NULL,
                    rangeRestr="either", expand=FALSE)
           standardGeneric("getTranscripts"))

setGeneric("getExons", signature=c("ann", "ranges"),
           function(ann, ranges=NULL, chromosome=NULL, strand=NULL,
                    rangeRestr="either", expand=FALSE)
           standardGeneric("getExons"))
