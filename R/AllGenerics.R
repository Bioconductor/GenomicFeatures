### Should we rename the first arg. 'txdb' or just 'x'?
### Also, should we change the signature to dispatch on the first arg. only?

setGeneric("getTranscripts", signature="ranges",
           function(ranges=NULL, ann, rangeRestr="either", expand=FALSE)
           standardGeneric("getTranscripts"))

setGeneric("mapTranscripts", signature="ranges",
           function(ranges=NULL, ann)
           standardGeneric("mapTranscripts"))

setGeneric("getExons", signature="ranges",
           function(ranges=NULL, ann, rangeRestr="either", expand=FALSE)
           standardGeneric("getExons"))

setGeneric("mapExons", signature="ranges",
           function(ranges=NULL, ann)
           standardGeneric("mapExons"))
