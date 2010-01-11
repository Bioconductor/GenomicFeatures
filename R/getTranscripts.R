

## Get the list of possible values for a type from a table:
.getValsFromTable <- function(txann, field, table){
  sql <- paste("SELECT",field,"FROM",table)
  as.character(unlist(unique(dbGetQuery(txann@conn, sql))))
}

## Check that the type of thing (chromosome, strand etc.) is in the
## transcripts table
.checkFields <- function(txann, type = c("chromosome","strand"), data, table){
  type <- match.arg(type)
  annot <- .getValsFromTable(txann, type, table) 
  all(data %in% annot)
}



## Helper function to construct the tail end of the queries
.makeRangeSQL <- function(start, end, ranges, rangeRestr) {
  switch(rangeRestr,
         "both"   = paste("(", start, " >= ", start(ranges),
                          " AND ", end, " <= ", end(ranges), ")",
                          sep = ""),
         "either" = paste("(", start, " >= ", start(ranges),
                          " OR ", end, " <= ", end(ranges), ")",
                          sep = ""),
         "start"  = paste("(", start," >= ", start(ranges), ")",
                          sep = ""),
         "end"    = paste("(", end, " <= ", end(ranges), ")",
                          sep = ""))
}




## Method for RangedData objects
setMethod("getTranscripts", "RangedData",
    function(ranges=NULL, ann, rangeRestr="either", expand=FALSE)
    {
      if(.checkFields(ann, "strand", ranges[["strand"]], "transcripts")){
        strand <- ranges[["strand"]]
      }else{stop("Strand values for ranges do not match annotation DB")}

      ## check that the chromosomes are what we expect. 
      if(.checkFields(ann, "chromosome", space(ranges), "transcripts")){
        chromosome <- space(ranges)
      }else{stop("Space values for ranges do not match annotation DB")}

      ranges <- unlist(ranges(ranges), use.names=FALSE)
      .getTranscripts(txdb=ann, ranges=ranges,
                      chromosome=chromosome, strand=strand,
                      rangeRestr=rangeRestr, expand=expand)
    }
)


## If there is not ranged Data object, then we just want it all...
setMethod("getTranscripts", "missing",
    function(ranges=NULL, ann, rangeRestr="either", expand=FALSE)
    {
      .getTranscripts(txdb=ann, rangeRestr=rangeRestr, expand=expand)
    }
)


## convenience function for seeing the SQL
.printSQL <- function(sql) {
  cat(strwrap(gsub("\\n +"," ",sql)),sep="\n")
}

### Extract selected transcripts from 'txdb'.
.getTranscripts <- function(txdb, ranges=NULL, chromosome=NULL,
                            strand=NULL,
                            rangeRestr=c("both","either","start","end"),
                            expand=FALSE) {
  rangeRestr <- match.arg(rangeRestr)
  len <- max(length(chromosome), length(strand), length(ranges))
  sql <- paste("SELECT t.tx_id, t.chromosome, t.strand,",
               "tt.tx_start, tt.tx_end",
               "FROM transcripts AS t,",
               "transcripts_rtree AS tt",
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
  ans <- dbGetQuery(txdb@conn, sql)
  ans <- RangedData(ranges     = IRanges(start = ans[["tx_start"]],
                                         end   = ans[["tx_end"]]),
                    strand     = strand(ans[["strand"]]),
                    transcript = ans[["tx_id"]],
                    space      = ans[["chromosome"]])
  if (expand) {
    if (len == 0) {
      sqlexons <- paste("SELECT ests._tx_id, et.exon_start, et.exon_end",
                        "FROM exons_rtree AS et,",
                        "exons_transcripts AS ests",
                        "WHERE et._exon_id=ests._exon_id")
    } else {
      sqlexons <- paste("SELECT ests._tx_id, et.exon_start, et.exon_end",
                        "FROM exons_rtree AS et,",
                        "exons_transcripts AS ests,",
                        "transcripts AS t,",
                        "transcripts_rtree AS tt",
                        "WHERE et._exon_id=ests._exon_id",
                        "AND ests._tx_id=t._tx_id",
                        "AND t._tx_id=tt._tx_id",
                        sqlwhere)
    }
    sqlexons <- paste(sqlexons, "ORDER BY ests._tx_id, ests.exon_rank")
    if (getOption("verbose", FALSE)) {
      .printSQL(sqlexons)
    }
    exons <- dbGetQuery(txdb@conn, sqlexons)
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




setMethod("getExons", "RangedData",
    function(ranges=NULL, ann, rangeRestr="either", expand=FALSE)
    {
      if(.checkFields(ann, "strand", ranges[["strand"]], "exons")){
        strand <- ranges[["strand"]]
      }else{stop("Strand values for ranges do not match annotation DB")}

      ## check that the chromosomes are what we expect. 
      if(.checkFields(ann, "chromosome", space(ranges), "exons")){
        chromosome <- space(ranges)
      }else{stop("Space values for ranges do not match annotation DB")}
      
      ranges <- unlist(ranges(ranges), use.names=FALSE)
      .getExons(txdb=ann, ranges=ranges,
                chromosome=chromosome, strand=strand,
                rangeRestr=rangeRestr, expand=expand)
    }
)

## No ranged Data object = get everything.
setMethod("getExons", "missing",
    function(ranges=NULL, ann, rangeRestr="either", expand=FALSE)
    {
      .getExons(txdb=ann, rangeRestr=rangeRestr, expand=expand)
    }
)



### Extract selected exons from 'txdb'.
.getExons <- function(txdb, ranges=NULL, chromosome=NULL,
                      strand=NULL,
                      rangeRestr=c("both","either","start","end"),
                      expand=FALSE) {
  rangeRestr <- match.arg(rangeRestr)
  len <- max(length(chromosome), length(strand), length(ranges))
  sql <- paste("SELECT e.chromosome, e.strand,",
               "et.exon_start, et.exon_end",
               "FROM exons AS e,",
               "exons_rtree AS et",
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
  ans <- dbGetQuery(txdb@conn, sql)
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
                     "exons_rtree AS et",
                     "WHERE ests._tx_id=t._tx_id",
                     "AND ests._exon_id=e._exon_id",
                     "AND e._exon_id=et._exon_id",
                     sqlwhere)
    }
    sqltx <- paste(sqltx, "ORDER BY ests._exon_id")
    if (getOption("verbose", FALSE)) {
      .printSQL(sqltx)
    }
    tx <- dbGetQuery(txdb@conn, sqltx)
    ans[["transcripts"]] <-
      IRanges:::newCompressedList("CompressedCharacterList",
                                  as.character(tx[["tx_id"]]),
                                  end = end(Rle(tx[["_exon_id"]])))
  }
  ans
}



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

