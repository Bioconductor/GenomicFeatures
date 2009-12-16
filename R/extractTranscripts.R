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

## put myTest.sqlite into /data and formalize the tests into unit tests.


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



##method to scan transcripts
getTranscripts <- function(transcript, ranges=NULL, chromosome=NULL,
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
##   sqlRange <- paste(" (tx_start>", start(ranges),
##                 " AND tx_end<", end(ranges),")")
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
  if(showSQL==TRUE){print(sql)}
  dbGetQuery(con, sql)
}


##TODO: for unit tests, put myTest.sqlite into /data and load it as needed.

## ##eg
## ## tx = loadFeatures("myTest.sqlite")

## ##This works though:
## ##foo = IRanges(start=c(500), end=c(10000))
## foo = IRanges(start=c(1000), end=c(20000))
## getTranscripts(tx,foo,"chr1","-")


## ##BUT NOT this (because no data there):
## foo = IRanges(start=c(16000), end=c(20000))
## getTranscripts(tx,foo,"chr1","-")


## ##AND a more complicated example...
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getTranscripts(tx,foo,"chr1","-")

##Compound search:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getTranscripts(tx,foo,c("chr1","chr2"),c("-","+"))

##Compound expanded search:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getTranscripts(tx,foo,c("chr1","chr2"),c("-","+"), expand=TRUE)



## ##method to scan exons
getExons <- function(transcript, ranges=NULL, chromosome=NULL,
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
##   sqlRange <- paste(" (exon_start>", start(ranges),
##                 " AND exon_end<", end(ranges),")")
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
  if(showSQL==TRUE){print(sql)}
  dbGetQuery(con, sql)  
}


## ##Complex example for exons.
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(tx,foo,"chr1","-")



##vectorized tests:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(tx,foo,c("chr1","chr2"),c("-","+"))

## expanded:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(tx,foo,c("chr1","chr2"),c("-","+"), expand=TRUE, showSQL=TRUE)
