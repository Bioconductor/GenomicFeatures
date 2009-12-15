## TODO:
## Chromosome and strand need to be vectorized like ranges. - DONE

## Append Chromosome, strand and ranges where clauses separately depending ... - DONE

## Add parameter to allow expanded query... (basically also join in exons when
## dealing with transcripts or vice versa for getExons, where we would have to
## also join in transcripts).

## expand to be a set of methods.  So this core will then become a helper
## function, and the methods will allow me to pass in either a RangedData
## object or a Ranges List.  Or other ways of representing the same kind of data.


##method to get transcripts
getTranscripts <- function(transcript, ranges=NULL, chromosome=NULL, strand=NULL){
  con = transcript@con
  sqlbase <- paste("SELECT t._tx_id, tx_id, chromosome,
                    strand, tx_start, tx_end
                    FROM transcripts AS t,
                    transcript_tree AS tt
                    WHERE t._tx_id=tt._tx_id", sep="")

  
  sqlChromosome <- paste(" t.chromosome='",chromosome,"' ",sep="")

  sqlStrand <- paste(" t.strand='",strand,"' ",sep="")
  
  sqlRange <- paste(" (tx_start>", start(ranges),
                " AND tx_end<", end(ranges),")")

  sql <- paste(sqlbase,
               " AND (",
               paste(sqlChromosome, collapse=" OR "),
               ") AND (",
               paste(sqlStrand, collapse=" OR "),
               ") AND (",
               paste(sqlRange,collapse=" OR "),
               ")")
  ##print(sql)
  dbGetQuery(con, sql)
}






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



## ##method to get exons
getExons <- function(transcript, ranges, chromosome, strand){
  con = transcript@con
  sqlbase <- paste("SELECT e._exon_id, chromosome, strand,
                    exon_start, exon_end
                    FROM exons AS e,
                    exon_tree AS et
                    WHERE e._exon_id=et._exon_id", sep="")

  sqlChromosome <- paste(" e.chromosome='",chromosome,"' ",sep="")

  sqlStrand <- paste(" e.strand='",strand,"' ",sep="")
  
  sqlRange <- paste(" (exon_start>", start(ranges),
                " AND exon_end<", end(ranges),")")

  sql <- paste(sqlbase,
               " AND (",
               paste(sqlChromosome, collapse=" OR "),
               ") AND (",
               paste(sqlStrand, collapse=" OR "),
               ") AND (",
               paste(sqlRange,collapse=" OR "),
               ")")
  ##print(sql)
  dbGetQuery(con, sql)  
}

## ##Complex example for exons.
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(tx,foo,"chr1","-")



##vectorized tests:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(tx,foo,c("chr1","chr2"),c("-","+"))
