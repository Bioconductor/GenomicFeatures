### Maybe we'll support other types of GenomicFeatures objects later.
setClass("GenomicFeatures", contains="VIRTUAL")

### Concrete GenomicFeatures type.
setClass("TranscriptDb",
    contains="GenomicFeatures",
    representation(envir="environment")
)

