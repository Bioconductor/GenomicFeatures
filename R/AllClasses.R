### Maybe we'll support other types of GenomicFeatures objects later.
#setClass("GenomicFeatures", contains="VIRTUAL")

############################################################################
## We define a new class when we have a new kind of schema to represent.

### Concrete GenomicFeatures types
.TranscriptDb <-
    setRefClass("TranscriptDb", contains="AnnotationDb",
        fields=list(isActiveSeq="logical"),
                methods=list(
                  initialize=function(){
                    .conn <-
                      if (missing(sqliteFile)) dbConnect(SQLite())
                      else dbConnect(SQLite(), sqliteFile)
                    ## then get default values for ActiveSeqs
                    seqNames <- .getChromInfo(conn)$chrom
                    seqNVals <- rep(TRUE, length(seqNames))
                    names(seqNVals) <- seqNames
                    callSuper(..., conn=.conn, isActiveSeq=seqNVals)
                  }))

.FeatureDb <-
    setRefClass("FeatureDb", contains="AnnotationDb")

