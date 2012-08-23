## transcripts.R
setGeneric("transcripts", function(x, ...) standardGeneric("transcripts"))

setGeneric("exons", function(x, ...) standardGeneric("exons"))

setGeneric("cds", function(x, ...) standardGeneric("cds"))


setGeneric("microRNAs", function(x) standardGeneric("microRNAs"))

setGeneric("tRNAs", function(x) standardGeneric("tRNAs"))

setGeneric("promoters", signature="x",
    function(x, upstream=2000, downstream=200, ...)
        standardGeneric("promoters"))

setGeneric("transcriptsBy", signature="x",
    function(x, by=c("gene", "exon", "cds"), ...)
        standardGeneric("transcriptsBy")
)

## transcriptsBy.R
setGeneric("exonsBy", signature="x",
    function(x, by=c("tx", "gene"), ...) standardGeneric("exonsBy")
)

setGeneric("cdsBy", signature="x",
    function(x, by=c("tx", "gene"), ...) standardGeneric("cdsBy")
)

setGeneric("intronsByTranscript",
    function(x, ...) standardGeneric("intronsByTranscript")
)

setGeneric("fiveUTRsByTranscript", 
    function(x, ...) standardGeneric("fiveUTRsByTranscript")
)

setGeneric("threeUTRsByTranscript", 
    function(x, ...) standardGeneric("threeUTRsByTranscript")
)

setGeneric("transcriptsByOverlaps", signature="x",
    function(x, ranges, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"), ...)
        standardGeneric("transcriptsByOverlaps")
)

## transciptsByOverlaps.R
setGeneric("exonsByOverlaps", signature="x",
    function(x, ranges, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"), ...)
        standardGeneric("exonsByOverlaps")
)

setGeneric("cdsByOverlaps", signature="x",
    function(x, ranges, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"), ...)
        standardGeneric("cdsByOverlaps")
)

## TranscriptDb-class.R
setGeneric("saveFeatures", signature="x",
           function(x, file) standardGeneric("saveFeatures"))

## features.R
setGeneric("features", signature="x",
           function(x) standardGeneric("features"))


## isActiveSeq
setGeneric("isActiveSeq", function(x) standardGeneric("isActiveSeq"))
setGeneric("isActiveSeq<-",function(x, value) standardGeneric("isActiveSeq<-"))


## generics for seqnamStyle and friends
setGeneric("determineDefaultSeqnameStyle",
           function(x) standardGeneric("determineDefaultSeqnameStyle"))
setGeneric("seqnameStyle<-", signature="x",
           function(x, value) standardGeneric("seqnameStyle<-"))

