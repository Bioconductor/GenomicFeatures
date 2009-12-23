### Plan to support other types of GenomicAnnotation objects later.
### HP: I suggest to either rename this class GenomicFeatures (or
### Features), or to rename saveFeatures()/loadFeatures() ->
### saveGenomicAnnotation()/loadGenomicAnnotation().
setClass("GenomicAnnotation",
         contains="VIRTUAL")

### Concrete GenomicAnnotation type.
setClass("TranscriptAnnotation",
         contains="GenomicAnnotation",
         representation=representation(conn = "SQLiteConnection"))

