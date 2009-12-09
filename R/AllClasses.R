## Plan to make other different Objects later.
setClass("GenomicAnnotation",
         contains="VIRTUAL")


## We need to define an actual Transcript Object.
setClass("TranscriptAnnotation",
         contains="GenomicAnnotation",
         representation=representation(con = "SQLiteConnection"))



