### Maybe we'll support other types of GenomicFeatures objects later.
#setClass("GenomicFeatures", contains="VIRTUAL")

############################################################################
## We define a new class when we have a new kind of schema to represent.

### Concrete GenomicFeatures types
setClass("TranscriptDb",
    contains="AnnotationDb",
    representation(envir="environment", isActiveSeq="logical")
)

setClass("FeatureDb",
    contains="AnnotationDb",
    representation(envir="environment")
)
