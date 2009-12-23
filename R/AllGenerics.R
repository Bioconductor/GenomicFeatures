setGeneric("getTranscripts", signature=c("ann", "ranges"),
           function(ann, ranges=NULL, chromosome=NULL, strand=NULL,
                    rangeRestr="either", expand=FALSE, showSQL=FALSE)
           standardGeneric("getTranscripts"))

setGeneric("getExons", signature=c("ann", "ranges"),
           function(ann, ranges=NULL, chromosome=NULL, strand=NULL,
                    rangeRestr="either", expand=FALSE, showSQL=FALSE)
           standardGeneric("getExons"))
