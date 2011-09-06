### Maybe we'll support other types of GenomicFeatures objects later.
#setClass("GenomicFeatures", contains="VIRTUAL")

############################################################################
## We define a new class when we have a new kind of schema to represent.

### Concrete GenomicFeatures types
.TranscriptDb <-
    setRefClass("TranscriptDb", contains="AnnotationDb",
        fields=list(isActiveSeq="logical"))

.FeatureDb <-
    setRefClass("FeatureDb", contains="AnnotationDb")
