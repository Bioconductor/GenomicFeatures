setGeneric("getTranscripts", signature="x",
           function(x, transcript, chromosome=NULL, strand=NULL,
                    rangeRestr="either", expand=FALSE, showSQL=FALSE)
           standardGeneric("getTranscripts"))

setGeneric("getExons", signature="x",
           function(x, transcript, chromosome=NULL, strand=NULL,
                    rangeRestr="either", expand=FALSE, showSQL=FALSE)
           standardGeneric("getExons"))
