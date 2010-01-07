### Should we rename the first arg. 'txdb' or just 'x'?
### Also, should we change the signature to dispatch on the first arg. only?

setGeneric("getTranscripts", signature=c("ann", "ranges"),
           function(ann, ranges=NULL, chromosome=NULL, strand=NULL,
                    rangeRestr="either", expand=FALSE)
           standardGeneric("getTranscripts"))

setGeneric("mapTranscripts", signature=c("ranges"),
           function(ranges=NULL, ann)
           standardGeneric("mapTranscripts"))

setGeneric("getExons", signature=c("ann", "ranges"),
           function(ann, ranges=NULL, chromosome=NULL, strand=NULL,
                    rangeRestr="either", expand=FALSE)
           standardGeneric("getExons"))
