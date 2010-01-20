

setMethod("transcriptsByRanges", "RangedData",
    function(txdb, ranges, restrict = "any", columns=c("tx_id", "tx_name"))
    {
      ## check that txdb is in fact a TranscriptDb object
      if(!is(txdb,"TranscriptDb"))stop("txdb MUST be a TranscriptDb object.")
      ## check that any supplied strands are what we expect.
      ## note that NULL is ok and will not trip this error.
      ## If there are NULL strands or chromosomes, then .mapTranscripts() will
      ## have to search a bit differently
      if(.checkFields(txdb, "tx_strand", ranges[["strand"]], "transcript")){
        strand <- ranges[["strand"]]
      }else{stop("Strand values for ranges do not match annotation DB")}
      ## check that the chromosomes are what we expect. 
      if(.checkFields(txdb, "tx_chrom", space(ranges), "transcript")){
        chrom <- space(ranges)
      }else{stop("Space values for ranges do not match annotation DB")}
      ranges <- unlist(ranges(ranges), use.names=FALSE)
      .mapTranscripts(txdb=txdb, ranges=ranges,
                      chrom=chrom, strand=strand,
                      type = restrict, col=columns,
                      format="get")
    }
)



## If there is not ranged Data object, then we just want it all...
setMethod("transcriptsByRanges", "missing",
    function(txdb, ranges, restrict = "any", columns=c("tx_id", "tx_name"))
    {
      .mapTranscripts(txdb=txdb, ranges=NULL,
                      chrom=NULL, strand=NULL,
                      type=restrict, col=columns,
                      format="get")
    }
)


##TODO: for unit tests, put myTest.sqlite into /data and load it as needed.

## ##eg
## ## library(GenomicFeatures)
## ## txdb = loadFeatures("myTest.sqlite")
## txdb <- loadFeatures(system.file("extdata", "HG18test.sqlite",
##                       package="GenomicFeatures"))


## ##This works though:
## ##foo = IRanges(start=c(500), end=c(10000))
## foo = IRanges(start=c(1000), end=c(20000))
## getTranscripts(txdb,foo,"chr1","-")


## ##BUT NOT this (because no data there):
## foo = IRanges(start=c(16000), end=c(20000))
## getTranscripts(txdb,foo,"chr1","-", rangeRestr="both")


## ##AND a more complicated example...
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getTranscripts(txdb,foo,"chr1","-")

##Compound search:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getTranscripts(txdb,foo,c("chr1","chr2"),c("-","+"))

##Compound expanded search:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getTranscripts(txdb,foo,c("chr1","chr2"),c("-","+"),
## rangeRestr= "both", expand=TRUE)


## Example for RangedData and chrom:
## rd1 <- RangedData(foo, space="chr1")
## rd2 <- RangedData(foo, space=c("chr2","chr3"))



setMethod("exonsByRanges", "RangedData",
    function(txdb, ranges, restrict = "any", columns=c("exon_id", "exon_name"))
    {
      ## check that txdb is in fact a TranscriptDb object
      if(!is(txdb,"TranscriptDb"))stop("txdb MUST be a TranscriptDb object.")
      ## check that any supplied strands are what we expect.
      ## note that NULL is ok and will not trip this error.
      ## If there are NULL strands or chromosomes, then .mapTranscripts() will
      ## have to search a bit differently
      if(.checkFields(txdb, "exon_strand", ranges[["exon_strand"]], "exon")){
        strand <- ranges[["exon_strand"]]
      }else{stop("Strand values for ranges do not match annotation DB")}
      ## check that the chromosomes are what we expect. 
      if(.checkFields(txdb, "exon_chrom", space(ranges), "exon")){
        chrom <- space(ranges)
      }else{stop("Space values for ranges do not match annotation DB")}
      ranges <- unlist(ranges(ranges), use.names=FALSE)
      .mapExons(txdb=txdb, ranges=ranges,
                chrom=chrom, strand=strand,
                type = restrict, col=columns,
                format="get")
    }
)


## If there is not ranged Data object, then we just want it all...
setMethod("exonsByRanges", "missing",
    function(txdb, ranges, restrict = "any", columns=c("exon_id", "exon_name"))
    {
      .mapExons(txdb=txdb, ranges=NULL,
                chrom=NULL, strand=NULL,
                type=restrict, col=columns,
                format="get")
    }
)




## ##Complex example for exons.
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(txdb,foo,"chr1","-")

##vectorized tests:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(txdb,foo,c("chr1","chr2"),c("-","+"))

## expanded:
## option(verbose = TRUE)
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(txdb,foo,c("chr1","chr2"),c("-","+"), expand=TRUE)


## LATER STUFF
## TODO:
## I need to add some discovery methods so that the users know what the values
## are that can be passed into chromosome etc. are.

## I will also need to add some warnings so that when users try to pass in a
## bad value for chromosome, they can get a message that says "your paramater
## needs to have values that look like: "chr1" "chr2" etc."


## Add the following test case (and modify the above test cases and add one or ttwo of those as well)

## library(GenomicFeatures)
## txdb <- loadFeatures("testDB.sqlite")
## library(IRanges)
## ranges = IRanges(start=c(500,10500), end=c(10000,30000))
## rd = RangedData(ranges=ranges, space = c("chr1","chr2"), c("-","-"))
## getExons(rd,txdb)

