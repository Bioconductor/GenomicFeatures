### Plan to support other types of GenomicAnnotation objects later.
setClass("GenomicAnnotation",
         contains="VIRTUAL")

### Concrete GenomicAnnotation type.
setClass("TranscriptAnnotation",
         contains="GenomicAnnotation",
         representation=representation(conn = "SQLiteConnection"))

