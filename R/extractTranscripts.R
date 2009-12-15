## TODO:
## Chromosome and strand need to be vectorized like ranges.

## Append Chromosome, strand and ranges where clauses separately depending ...

## Add parameter to allow expanded query... (basically also join in exons when
## dealing with transcripts or vice versa for getExons, where we would have to
## also join in transcripts).

## expand to be a set of methods.  So this core will then become a helper
## function, and the methods will allow me to pass in either a RangedData
## object or a Ranges List.  Or other ways of representing the same kind of data.


##method to get transcripts
getTranscripts <- function(transcript, ranges=NULL, chromosome=NULL, strand=NULL){
  ## validObject(ranges) ##How to test? (can't seem to make a "bad" one)
  con = transcript@con
  sql1 <- paste("SELECT t._tx_id, tx_id, chromosome,
                strand, tx_start, tx_end
                FROM transcripts AS t,
                transcript_tree AS tt
                WHERE t._tx_id=tt._tx_id
                AND t.chromosome='",chromosome,
                "' AND t.strand='",strand,"' AND (",sep="")

  sql2 <- paste(" (tx_start>", start(ranges),
                " AND tx_end<", end(ranges),")")

  sql <- paste(sql1, paste(sql2,collapse=" OR "), " )")
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




## ##method to get exons
getExons <- function(transcript, ranges, chromosome, strand){
  ## validObject(ranges) ##How to test?
  con = transcript@con
  sql1 <- paste("SELECT e._exon_id, chromosome, strand,
                exon_start, exon_end
                FROM exons AS e,
                exon_tree AS ee
                WHERE e._exon_id=ee._exon_id
                AND e.chromosome='",chromosome,
                "' AND e.strand='",strand,"' AND (",sep="")

  sql2 <- paste(" (exon_start>", start(ranges),
                " AND exon_end<", end(ranges),")")

  sql <- paste(sql1, paste(sql2,collapse=" OR "), " )")
  ##print(sql)
  dbGetQuery(con, sql)
}

## ##Complex example for exons.
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(tx,foo,"chr1","-")
