### Plan to support other types of GenomicFeatures objects later.
setClass("GenomicFeatures",
         contains="VIRTUAL")

### Concrete GenomicFeatures type.
setClass("TranscriptDb",
         contains="GenomicFeatures",
         representation=representation(conn = "SQLiteConnection"))

