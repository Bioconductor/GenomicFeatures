### Maybe we'll support other types of GenomicFeatures objects later.
setClass("GenomicFeatures", contains="VIRTUAL")

### Concrete GenomicFeatures types
setClass("TranscriptDb",
    contains="GenomicFeatures",
    representation(envir="environment")
)

setClass("FeatureDb",
    contains="GenomicFeatures",
    representation(envir="environment")
)
