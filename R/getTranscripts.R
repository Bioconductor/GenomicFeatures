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
## data. - DONE

## make a "mini" myTest.sqlite and then put it in... - DONE

## Now I can do this:
##tx = loadFeatures(system.file("inst/extdata/HG18.sqlite",package="GenomicFeatures"))

## put myTest.sqlite into /data and formalize the tests into unit tests.


## Helper function to construct the tail end of the queries
.makeRangeSQL <- function(start, end, ranges, rangeRestr) {
  switch(rangeRestr,
         "both"   = paste("(", start, " > ", start(ranges),
                          " AND ", end, " < ", end(ranges), ")",
                          sep = ""),
         "either" = paste("(", start, " > ", start(ranges),
                          " OR ", end, " < ", end(ranges), ")",
                          sep = ""),
         "start"  = paste("(", start," > ", start(ranges), ")",
                          sep = ""),
         "end"    = paste("(", end, " < ", end(ranges), ")",
                          sep = ""))
}


## getTranscript methods accomodate different object and try to use the
## information that they contain
setMethod("getTranscripts", c("TranscriptAnnotation", "IRanges"),
    function(ann, ranges=NULL, chromosome=NULL, strand=NULL,
             rangeRestr="either", expand=FALSE)
    {
      .getTranscripts(txann=ann, ranges=ranges,
                      chromosome=chromosome, strand=strand,
                      rangeRestr=rangeRestr, expand=expand)
    }
)


## If they give me a RangedData object or a RangesList object, then I should
## just check for any additional chromosome and strand information in there,
## and append that to my chromosome and strand information.

setMethod("getTranscripts", c("TranscriptAnnotation", "RangedData"),
    function(ann, ranges=NULL, chromosome=NULL, strand=NULL,
             rangeRestr="either", expand=FALSE)
    {
      rdChromosome <- space(ranges)
      chromosome <- c(chromosome, rdChromosome)
      ranges <- unlist(ranges(ranges), use.names=FALSE)
      .getTranscripts(txann=ann, ranges=ranges,
                      chromosome=chromosome, strand=strand,
                      rangeRestr=rangeRestr, expand=expand)
    }
)

setMethod("getTranscripts", c("TranscriptAnnotation", "RangesList"),
    function(ann, ranges=NULL, chromosome=NULL, strand=NULL,
             rangeRestr="either", expand=FALSE)
    {
      ranges <- unlist(ranges, use.names=FALSE)
      .getTranscripts(txann=ann, ranges=ranges,
                      chromosome=chromosome, strand=strand,
                      rangeRestr=rangeRestr, expand=expand)
    }
)


.printSQL <- function(sql) {
  cat(strwrap(gsub("\\n +"," ",sql)),sep="\n")
}

### Extract selected transcripts from 'txann'.
.getTranscripts <- function(txann, ranges=NULL, chromosome=NULL,
                            strand=NULL,
                            rangeRestr=c("both","either","start","end"),
                            expand=FALSE) {
  rangeRestr <- match.arg(rangeRestr)
  len <- max(length(chromosome), length(strand), length(ranges))
  sql <- paste("SELECT t.tx_id, t.chromosome, t.strand,",
               "tt.tx_start, tt.tx_end",
               "FROM transcripts AS t,",
               "transcript_tree AS tt",
               "WHERE t._tx_id=tt._tx_id")
  if (len > 0) {
    if (!is.null(chromosome)) {
      if (length(chromosome) < len)
        chromosome <- rep(chromosome, length.out = len)
      sqlwhere <- paste("t.chromosome='", chromosome, "'", sep="")
    } else {
      sqlwhere <- character(0)
    }
    if (!is.null(strand)) {
      if (length(strand) < len)
        strand <- rep(strand, length.out = len)
      sqladd <- paste("t.strand='", strand, "'", sep="")
      if (length(sqlwhere) == 0)
        sqlwhere <- sqladd
      else
        sqlwhere <- paste(sqlwhere, sqladd, sep = " AND ")
    }
    if (!is.null(ranges)) {
      if (length(ranges) < len)
        ranges <- rep(ranges, length.out = len)
      sqladd <- .makeRangeSQL("tt.tx_start", "tt.tx_end", ranges, rangeRestr)
      if (length(sqlwhere) == 0)
        sqlwhere <- sqladd
      else
        sqlwhere <- paste(sqlwhere, sqladd, sep = " AND ")
    }
    sqlwhere <- paste("AND (",
                      paste("(", sqlwhere, ")", sep = "", collapse = " OR "),
                      ")", sep = "")
    sql <- paste(sql, sqlwhere)
  }
  sql <- paste(sql, "ORDER BY t._tx_id")
  if (getOption("verbose", FALSE)) {
    .printSQL(sql)
  }
  ans <- dbGetQuery(txann@conn, sql)
  ans <- RangedData(ranges     = IRanges(start = ans[["tx_start"]],
                                         end   = ans[["tx_end"]]),
                    strand     = strand(ans[["strand"]]),
                    transcript = ans[["tx_id"]],
                    space      = ans[["chromosome"]])
  if (expand) {
    if (len == 0) {
      sqlexons <- paste("SELECT ests._tx_id, et.exon_start, et.exon_end",
                        "FROM exon_tree AS et,",
                        "exons_transcripts AS ests",
                        "WHERE et._exon_id=ests._exon_id")
    } else {
      sqlexons <- paste("SELECT ests._tx_id, et.exon_start, et.exon_end",
                        "FROM exon_tree AS et,",
                        "exons_transcripts AS ests,",
                        "transcripts AS t,",
                        "transcript_tree AS tt",
                        "WHERE et._exon_id=ests._exon_id",
                        "AND ests._tx_id=t._tx_id",
                        "AND t._tx_id=tt._tx_id",
                        sqlwhere)
    }
    sqlexons <- paste(sqlexons, "ORDER BY ests._tx_id, ests.exon_rank")
    if (getOption("verbose", FALSE)) {
      .printSQL(sqlexons)
    }
    exons <- dbGetQuery(txann@conn, sqlexons)
    ans[["exons"]] <-
      IRanges:::newCompressedList("CompressedIRangesList",
                                  IRanges(start = exons[["exon_start"]],
                                          end   = exons[["exon_end"]]),
                                  end = end(Rle(exons[["_tx_id"]])))
  }
  ans
}


##TODO: for unit tests, put myTest.sqlite into /data and load it as needed.

## ##eg
## ## library(GenomicFeatures)
## ## txann = loadFeatures("myTest.sqlite")
## txann <- loadFeatures(system.file("extdata", "HG18test.sqlite",
##                       package="GenomicFeatures"))


## ##This works though:
## ##foo = IRanges(start=c(500), end=c(10000))
## foo = IRanges(start=c(1000), end=c(20000))
## getTranscripts(txann,foo,"chr1","-")


## ##BUT NOT this (because no data there):
## foo = IRanges(start=c(16000), end=c(20000))
## getTranscripts(txann,foo,"chr1","-", rangeRestr="both")


## ##AND a more complicated example...
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getTranscripts(txann,foo,"chr1","-")

##Compound search:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getTranscripts(txann,foo,c("chr1","chr2"),c("-","+"))

##Compound expanded search:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getTranscripts(txann,foo,c("chr1","chr2"),c("-","+"),
## rangeRestr= "both", expand=TRUE)


## Example for RangedData and chrom:
## rd1 <- RangedData(foo, space="chr1")
## rd2 <- RangedData(foo, space=c("chr2","chr3"))


setMethod("getExons", c("TranscriptAnnotation", "IRanges"),
    function(ann, ranges=NULL, chromosome=NULL, strand=NULL,
             rangeRestr="either", expand=FALSE)
    {
      .getExons(txann=ann, ranges=ranges,
                chromosome=chromosome, strand=strand,
                rangeRestr=rangeRestr, expand=expand)
    }
)

setMethod("getExons", c("TranscriptAnnotation", "RangedData"),
    function(ann, ranges=NULL, chromosome=NULL, strand=NULL,
             rangeRestr="either", expand=FALSE)
    {
      rdChromosome <- space(ranges)
      chromosome <- c(chromosome, rdChromosome)
      ranges <- unlist(ranges(ranges), use.names=FALSE)
      .getExons(txann=ann, ranges=ranges,
                chromosome=chromosome, strand=strand,
                rangeRestr=rangeRestr, expand=expand)
    }
)

setMethod("getExons", c("TranscriptAnnotation", "RangesList"),
    function(ann, ranges, chromosome=NULL,
             strand=NULL, rangeRestr="either",
             expand=FALSE)
    {
      ranges <- unlist(ranges, use.names=FALSE)
      .getExons(txann=ann, ranges=ranges,
                chromosome=chromosome, strand=strand,
                rangeRestr=rangeRestr, expand=expand)            
    }
)

### Extract selected exons from 'txann'.
.getExons <- function(txann, ranges=NULL, chromosome=NULL,
                      strand=NULL,
                      rangeRestr=c("both","either","start","end"),
                      expand=FALSE) {
  rangeRestr <- match.arg(rangeRestr)
  len <- max(length(chromosome), length(strand), length(ranges))
  sql <- paste("SELECT e.chromosome, e.strand,",
               "et.exon_start, et.exon_end",
               "FROM exons AS e,",
               "exon_tree AS et",
               "WHERE e._exon_id=et._exon_id")
  if (len > 0) {
    if (!is.null(chromosome)) {
      if (length(chromosome) < len)
        chromosome <- rep(chromosome, length.out = len)
      sqlwhere <- paste("e.chromosome='", chromosome, "'", sep="")
    } else {
      sqlwhere <- character(0)
    }
    if (!is.null(strand)) {
      if (length(strand) < len)
        strand <- rep(strand, length.out = len)
      sqladd <- paste("e.strand='", strand, "'", sep="")
      if (length(sqlwhere) == 0)
        sqlwhere <- sqladd
      else
        sqlwhere <- paste(sqlwhere, sqladd, sep = " AND ")
    }
    if (!is.null(ranges)) {
      if (length(ranges) < len)
        ranges <- rep(ranges, length.out = len)
      sqladd <- .makeRangeSQL("et.exon_start", "et.exon_end", ranges, rangeRestr)
      if (length(sqlwhere) == 0)
        sqlwhere <- sqladd
      else
        sqlwhere <- paste(sqlwhere, sqladd, sep = " AND ")
    }
    sqlwhere <- paste("AND (",
                      paste("(", sqlwhere, ")", sep = "", collapse = " OR "),
                      ")", sep = "")
    sql <- paste(sql, sqlwhere)
  }
  sql <- paste(sql, "ORDER BY e._exon_id")
  if (getOption("verbose", FALSE)) {
    .printSQL(sql)
  }
  ans <- dbGetQuery(txann@conn, sql)
  ans <- RangedData(ranges = IRanges(start = ans[["exon_start"]],
                                     end   = ans[["exon_end"]]),
                    strand = strand(ans[["strand"]]),
                    space  = ans[["chromosome"]])
  if (expand) {
    if (len == 0) {
      sqltx <- paste("SELECT ests._exon_id, t.tx_id",
                     "FROM exons_transcripts AS ests,",
                     "transcripts AS t",
                     "WHERE ests._tx_id=t._tx_id")
    } else {
      sqltx <- paste("SELECT ests._exon_id, t.tx_id",
                     "FROM exons_transcripts AS ests,",
                     "transcripts AS t,",
                     "exons AS e,",
                     "exon_tree AS et",
                     "WHERE ests._tx_id=t._tx_id",
                     "AND ests._exon_id=e._exon_id",
                     "AND e._exon_id=et._exon_id",
                     sqlwhere)
    }
    sqltx <- paste(sqltx, "ORDER BY ests._exon_id")
    if (getOption("verbose", FALSE)) {
      .printSQL(sqltx)
    }
    tx <- dbGetQuery(txann@conn, sqltx)
    ans[["transcripts"]] <-
      IRanges:::newCompressedList("CompressedCharacterList",
                                  as.character(tx[["tx_id"]]),
                                  end = end(Rle(tx[["_exon_id"]])))
  }
  ans
}


## ##Complex example for exons.
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(txann,foo,"chr1","-")

##vectorized tests:
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(txann,foo,c("chr1","chr2"),c("-","+"))

## expanded:
## option(verbose = TRUE)
## foo = IRanges(start=c(500,10500), end=c(10000,30000))
## getExons(txann,foo,c("chr1","chr2"),c("-","+"), expand=TRUE)


## LATER STUFF
## TODO:
## I need to add some discovery methods so that the users know what the values
## are that can be passed into chromosome etc. are.

## I will also need to add some warnings so that when users try to pass in a
## bad value for chromosome, they can get a message that says "your paramater
## needs to have values that look like: "chr1" "chr2" etc."
