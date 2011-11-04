### Maybe we'll support other types of GenomicFeatures objects later.
#setClass("GenomicFeatures", contains="VIRTUAL")

############################################################################
## We define a new class when we have a new kind of schema to represent.

## This is to try and tidy up before setRefClass()
gc()

### Concrete GenomicFeatures types
.TranscriptDb <-
    setRefClass("TranscriptDb", contains="AnnotationDb",
        fields=list(isActiveSeq="logical"),
        methods=list(
          initialize=function(...) {
              callSuper(...)
              if (0L == length(dbListTables(conn))) {
                  .self$isActiveSeq <- logical()
              } else {
                  seqNames <- .getChromInfo(conn)$chrom
                  .self$isActiveSeq <-
                      structure(!logical(length(seqNames)), .Names=seqNames)
              }
          .self
      }))

.FeatureDb <-
    setRefClass("FeatureDb", contains="AnnotationDb")

## This is to try and tidy up after setRefClass()
#gc()
