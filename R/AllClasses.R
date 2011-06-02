### Maybe we'll support other types of GenomicFeatures objects later.
setClass("GenomicFeatures", contains="VIRTUAL")

#####################################################################
## the following looks odd.  Why have more than one class?  The
## distinction is the DB schema held behind each class

### Concrete GenomicFeatures types
setClass("TranscriptDb",
    contains="GenomicFeatures",
    representation(envir="environment", activeSeqs="logical")
    ## representation(envir="environment")
)

setClass("FeatureDb",
    contains="GenomicFeatures",
    representation(envir="environment")
)
