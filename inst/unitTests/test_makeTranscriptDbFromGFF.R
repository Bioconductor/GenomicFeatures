## Test the helper functions, and then test the big one.
## Unlike the other test methods, for this one we don't need to be testing
## that the source is still the same since we point to a file type instead of
## to a service

## 1st set up some resources for all the tests
## gffFile <- system.file("extdata","a.gff3",package="GenomicFeatures")
## gffformat <- "gff3"
## gtfFile <- system.file("extdata","Aedes_aegypti.partial.gtf",
##                        package="GenomicFeatures")
## gffformat <- "gtf"



## test .deduceExonRankings



## test .mergeFramesViaRanges


## test .deduceTranscriptsFromGTF


## test  .prepareGTFTables


## test  .prepareGFF3Tables




## test makeTranscriptDbFromGFF
## run and compare to existing DB (saved)
## like this
    ## ## compare
    ## ok <- GenomicFeatures:::compareTranscriptDbs(txdb1, txdb0)
    ## checkTrue(ok)

