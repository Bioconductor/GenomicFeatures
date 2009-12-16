## TODO:
## Chromosome and strand need to be vectorized like ranges. - DONE

## Append Chromosome, strand and ranges where clauses separately depending
## ... - DONE

## Add parameters to allow expanded querys... (basically also join in exons when
## dealing with transcripts or vice versa for getExons, where we would have to
## also join in transcripts). - DONE

## Add checks in case someone wants to search "all Chromosomes" or all
## ranges - so that they can leave it blank to do this - DONE

## ALSO: add another parameter with a set of values () to control whether the
## ranges have to include the feature, or simply overlap with it.  range =
## c("start", "both", "end", "either").  Where "both" is the default and
## requires that both ends of the range are contained within the start and end
## elements of an annotation, "either" would say (foo >X OR bar <Y), meaning
## that either thing can be ok (for overlaps), and start and end would only
## check one thing or the other... - DONE





## expand to be a set of methods.  So this core will then become a helper
## function, and the methods will allow me to pass in either a RangedData
## object or a Ranges List.  Or other ways of representing the same kind of
## data.


## make a "mini" myTest.sqlite and then put it in...
## put myTest.sqlite into /data and formalize the tests into unit tests.




## Helper function to construct the tail end of the queries
.makeRangeSQL <- function(start, end, ranges, rangeRestr){
  if(rangeRestr == "both"){
    sqlRange <- paste(" (",start,">", start(ranges),
                      " AND ",end,"<", end(ranges),")")
  }else if(rangeRestr == "either"){
    sqlRange <- paste(" (",start,">", start(ranges),
                      " OR ",end,"<", end(ranges),")")

  }else if(rangeRestr == "start"){
    sqlRange <- paste(" (",start,">", start(ranges),")")

  }else if(rangeRestr == "end"){
    sqlRange <- paste(" (",end,"<", end(ranges),")")

  }
}



## getTranscript methods accomodate different object and try to use the
## information that they contain
setMethod("getTranscripts", "IRanges",
    function(x, transcript, chromosome=NULL,
             strand=NULL, rangeRestr="either",
             expand=FALSE, showSQL=FALSE)
    {
      .getTranscripts(transcript=transcript, ranges=x,
                      chromosome=chromosome, strand=strand,
                      rangeRestr=rangeRestr, expand=expand,
                      showSQL=showSQL)      
    }
)


## If they give me a RangedData object or a RangesList object, then I should
## just check for any additional chromosome and strand information in there,
## and append that to my chromosome and strand information.

setMethod("getTranscripts", "RangedData",
    function(x, transcript, chromosome=NULL,
             strand=NULL, rangeRestr="either",
             expand=FALSE, showSQL=FALSE)
    {
      rdChromosome=space(x) ##Q: is this REALLY Safe??? "space()"? Really???
      chromosome = c(chromosome, rdChromosome)
      ranges = unlist(ranges(x))
      .getTranscripts(transcript=transcript, ranges=ranges,
                      chromosome=chromosome, strand=strand,
                      rangeRestr=rangeRestr, expand=expand,
                      showSQL=showSQL)      
    }
)

setMethod("getTranscripts", "RangesList",
    function(x, transcript, chromosome=NULL,
             strand=NULL, rangeRestr="either",
             expand=FALSE, showSQL=FALSE)
    {
      ranges = unlist(ranges(x))
      .getTranscripts(transcript=transcript, ranges=ranges,
                      chromosome=chromosome, strand=strand,
                      rangeRestr=rangeRestr, expand=expand,
                      showSQL=showSQL)            
    }
)



##method to scan transcripts
.getTranscripts <- function(transcript, ranges=NULL, chromosome=NULL,
                           strand=NULL,
                           rangeRestr=c("both","either","start","end"),
                           expand=FALSE, showSQL=FALSE){
  con = transcript@con
  rangeRestr <- match.arg(rangeRestr)
  if(expand==FALSE){
    sqlbase <- "SELECT t._tx_id, tx_id, chromosome,
                strand, tx_start, tx_end
                FROM transcripts AS t,
                transcript_tree AS tt
                WHERE t._tx_id=tt._tx_id"
  }else{
    sqlbase <- "SELECT t._tx_id, tx_id, chromosome,
                strand, tx_start, tx_end, ests._exon_id,
                exon_start, exon_end
                FROM transcripts AS t,
                transcript_tree AS tt,
                exons_transcripts AS ests,
                exon_tree AS et
                WHERE t._tx_id=tt._tx_id
                AND t._tx_id=ests._tx_id
                AND ests._exon_id=et._exon_id"
  }
  if(!is.null(chromosome)){
  sqlChromosome <- paste(" t.chromosome='",chromosome,"' ",sep="")
  }else{sqlChromosome<-"t._tx_id=t._tx_id"}
  
  if(!is.null(strand)){
  sqlStrand <- paste(" t.strand='",strand,"' ",sep="")
  }else{sqlStrand<-"t._tx_id=t._tx_id"}
  
  if(!is.null(ranges)){
    sqlRange <- .makeRangeSQL("tx_start", "tx_end", ranges, rangeRestr)
  }else{sqlRange<-"t._tx_id=t._tx_id"}

  sql <- paste(sqlbase,
               " AND (",
               paste(sqlChromosome, collapse=" OR "),
               ") AND (",
               paste(sqlStrand, collapse=" OR "),
               ") AND (",
               paste(sqlRange,collapse=" OR "),
               ") ")
  if(showSQL==TRUE){print(gsub("\\n               ","",sql))}
  dbGetQuery(con, sql)
}


##TODO: for unit tests, put myTest.sqlite into /data and load it as needed.

## ##eg
## ## library(GenomicFeatures)
## ## tx = loadFeatures("myTest.sqlite")

## ##This works though:
## ##foo = IRanges(start=c(500), end=c(10000))
## foo = IRanges(start=c(1000), end=c(20000))
## getTranscripts(foo,tx,"chr1","-")


## ##BUT NOT this (because no data there):
## foo = IRanges(start=c(16000), end=c(20000))
## getTranscripts(foo,tx,"chr1","-", rangeRestr="both")


## ##AND a more complicated example...
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getTranscripts(foo,tx,"chr1","-")

##Compound search:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getTranscripts(foo,tx,c("chr1","chr2"),c("-","+"))

##Compound expanded search:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getTranscripts(foo,tx,c("chr1","chr2"),c("-","+"),
## rangeRestr= "both", expand=TRUE)


## Example for RangedData and chrom:
## rd1 <- RangedData(foo, space="chr1")
## rd2 <- RangedData(foo, space=c("chr2","chr3"))










setMethod("getExons", "IRanges",
    function(x, transcript, chromosome=NULL,
             strand=NULL, rangeRestr="either",
             expand=FALSE, showSQL=FALSE)
    {
      .getExons(transcript=transcript, ranges=x,
                      chromosome=chromosome, strand=strand,
                      rangeRestr=rangeRestr, expand=expand,
                      showSQL=showSQL)      
    }
)

setMethod("getExons", "RangedData",
    function(x, transcript, chromosome=NULL,
             strand=NULL, rangeRestr="either",
             expand=FALSE, showSQL=FALSE)
    {
      rdChromosome=space(x) ##Q: is this REALLY Safe??? "space()"? Really???
      chromosome = c(chromosome, rdChromosome)
      ranges = unlist(ranges(x))
      .getExons(transcript=transcript, ranges=ranges,
                      chromosome=chromosome, strand=strand,
                      rangeRestr=rangeRestr, expand=expand,
                      showSQL=showSQL)      
    }
)

setMethod("getExons", "RangesList",
    function(x, transcript, chromosome=NULL,
             strand=NULL, rangeRestr="either",
             expand=FALSE, showSQL=FALSE)
    {
      ranges = unlist(ranges(x))
      .getExons(transcript=transcript, ranges=ranges,
                      chromosome=chromosome, strand=strand,
                      rangeRestr=rangeRestr, expand=expand,
                      showSQL=showSQL)            
    }
)


## ##method to scan exons
.getExons <- function(transcript, ranges=NULL, chromosome=NULL,
                     strand=NULL,
                     rangeRestr=c("both","either","start","end"),
                     expand=FALSE, showSQL=FALSE){
  con = transcript@con
  rangeRestr <- match.arg(rangeRestr)
  if(expand==FALSE){
    sqlbase <- "SELECT e._exon_id, chromosome, strand,
                exon_start, exon_end
                FROM exons AS e,
                exon_tree AS et
                WHERE e._exon_id=et._exon_id"
  }else{
    sqlbase <- "SELECT e._exon_id, e.chromosome, e.strand,
                exon_start, exon_end, ests._tx_id,
                tx_id
                FROM exons AS e,
                exon_tree AS et,
                exons_transcripts AS ests,
                transcripts AS t
                WHERE e._exon_id=et._exon_id
                AND e._exon_id=ests._exon_id
                AND ests._tx_id=t._tx_id"
  } 

  if(!is.null(chromosome)){
  sqlChromosome <- paste(" e.chromosome='",chromosome,"' ",sep="")
  }else{sqlChromosome<-"e._exon_id=e._exon_id"}

  if(!is.null(strand)){
  sqlStrand <- paste(" e.strand='",strand,"' ",sep="")
  }else{sqlStrand<-"e._exon_id=e._exon_id"}
  
  if(!is.null(ranges)){
    sqlRange <- .makeRangeSQL("exon_start", "exon_end", ranges, rangeRestr)  
  }else{sqlRange<-"e._exon_id=e._exon_id"}

  sql <- paste(sqlbase,
               " AND (",
               paste(sqlChromosome, collapse=" OR "),
               ") AND (",
               paste(sqlStrand, collapse=" OR "),
               ") AND (",
               paste(sqlRange,collapse=" OR "),
               ")")
  if(showSQL==TRUE){print(gsub("\\n               ","",sql))}
  dbGetQuery(con, sql)  
}


## ##Complex example for exons.
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(foo,tx,"chr1","-")

##vectorized tests:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(foo,tx,c("chr1","chr2"),c("-","+"))

## expanded:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(foo,tx,c("chr1","chr2"),c("-","+"), expand=TRUE, showSQL=TRUE)









## LATER STUFF
## TODO:
## I need to add some discovery methods so that the users know what the values
## are that can be passed into chromosome etc. are.

## I will also need to add some warnings so that when users try to pass in a
## bad value for chromosome, they can get a message that says "your paramater
## needs to have values that look like: "chr1" "chr2" etc."
