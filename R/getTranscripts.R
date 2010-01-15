

## Get the list of possible values for a type from a table:
.getValsFromTable <- function(txann, colname, table){
  sql <- paste("SELECT", colname, "FROM", table)
  as.character(unlist(unique(dbGetQuery(txann@conn, sql))))
}

## Check that the type of thing (chromosome, strand etc.) is in the
## transcript table
.checkFields <- function(txann, colname, data, table){
  annot <- .getValsFromTable(txann, colname, table) 
  all(data %in% annot)
}



## Helper function to construct the tail end of the queries
.makeRangeSQL <- function(start, end, ranges, rangeRestr) {
  switch(rangeRestr,
         "both"   = paste("(", start, " >= ", start(ranges),
                          " AND ", end, " <= ", end(ranges), ")",
                          sep = ""),
         "either" = paste("((", start, " >= ", start(ranges),
                          " AND ", start, " <= ", end(ranges),")",
                          " OR (", end, " >= ", start(ranges),
                          " AND ", end, " <= ", end(ranges), ")",
                          " OR (", start, " <= ", start(ranges),
                          " AND ", end, " >= ", end(ranges), "))",
                          sep = ""),
         "start"  = paste("(", start, " >= ", start(ranges),
                          " AND ", start, " <= ", end(ranges), ")",
                          sep = ""),
         "end"    = paste("(", end, " >= ", start(ranges),
                          " AND ", end, " <= ", end(ranges), ")",
                          sep = "")) 
}




## Method for RangedData objects
setMethod("getTranscripts", "RangedData",
    function(ranges=NULL, ann, rangeRestr="either", expand=FALSE)
    {
      if(.checkFields(ann, "tx_strand", ranges[["tx_strand"]], "transcript")){
        strand <- ranges[["tx_strand"]]
      }else{stop("Strand values for ranges do not match annotation DB")}

      ## check that the chromosomes are what we expect. 
      if(.checkFields(ann, "tx_chrom", space(ranges), "transcript")){
        chrom <- space(ranges)
      }else{stop("Space values for ranges do not match annotation DB")}

      ranges <- unlist(ranges(ranges), use.names=FALSE)
      .getTranscripts(txdb=ann, ranges=ranges,
                      chrom=chrom, strand=strand,
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
.getTranscripts <- function(txdb, ranges=NULL, chrom=NULL,
                            strand=NULL,
                            rangeRestr=c("both","either","start","end"),
                            expand=FALSE) {
  rangeRestr <- match.arg(rangeRestr)
  len <- max(length(chrom), length(strand), length(ranges))
  ## Note that .mapTranscripts() uses the same SQL query.
  sql <- paste("SELECT tx_id, tx_name, tx_chrom, tx_strand,",
               "tx_start, tx_end",
               "FROM transcript INNER JOIN transcript_rtree",
               "ON (transcript._tx_id=transcript_rtree._tx_id)")
  if (len > 0) {
    if (!is.null(chrom)) {
      if (length(chrom) < len)
        chrom <- rep(chrom, length.out = len)
      sqlwhere <- paste("tx_chrom='", chrom, "'", sep="")
    } else {
      sqlwhere <- character(0)
    }
    if (!is.null(strand)) {
      if (length(strand) < len)
        strand <- rep(strand, length.out = len)
      sqladd <- paste("tx_strand='", strand, "'", sep="")
      if (length(sqlwhere) == 0)
        sqlwhere <- sqladd
      else
        sqlwhere <- paste(sqlwhere, sqladd, sep = " AND ")
    }
    if (!is.null(ranges)) {
      if (length(ranges) < len)
        ranges <- rep(ranges, length.out = len)
      sqladd <- .makeRangeSQL("tx_start", "tx_end", ranges, rangeRestr)
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
  sql <- paste(sql, "ORDER BY tx_id")
  if (getOption("verbose", FALSE)) {
    .printSQL(sql)
  }
  ans <- dbGetQuery(txdb@conn, sql)
  ans <- RangedData(ranges     = IRanges(start = ans[["tx_start"]],
                                         end   = ans[["tx_end"]]),
                    strand     = strand(ans[["tx_strand"]]),
                    GF_txId    = ans[["tx_id"]],
                    txName = ans[["tx_name"]], ## temp. just for troubleshooting
                    space      = ans[["tx_chrom"]])
  if (expand) {
    sqlexons <- paste("SELECT tx_id, exon_start, exon_end",
                      "FROM transcript",
                      "INNER JOIN splicing",
                      "ON (transcript._tx_id=splicing._tx_id)",
                      "INNER JOIN exon_rtree",
                      "ON (splicing._exon_id=exon_rtree._exon_id)")
    if (len != 0)
      sqlexons <- paste(sqlexons,
                        "INNER JOIN transcript_rtree",
                        "ON (transcript._tx_id=transcript_rtree._tx_id)",
                        "WHERE", sqlwhere)
    sqlexons <- paste(sqlexons, "ORDER BY tx_id, exon_rank")
    if (getOption("verbose", FALSE)) {
      .printSQL(sqlexons)
    }
    exons <- dbGetQuery(txdb@conn, sqlexons)
    ans[["exon"]] <-
      IRanges:::newCompressedList("CompressedIRangesList",
                                  IRanges(start = exons[["exon_start"]],
                                          end   = exons[["exon_end"]]),
                                  end = end(Rle(exons[["tx_id"]])))
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
      if(.checkFields(ann, "exon_strand", ranges[["exon_strand"]], "exon")){
        strand <- ranges[["exon_strand"]]
      }else{stop("Strand values for ranges do not match annotation DB")}

      ## check that the chromosomes are what we expect. 
      if(.checkFields(ann, "exon_chrom", space(ranges), "exon")){
        chrom <- space(ranges)
      }else{stop("Space values for ranges do not match annotation DB")}
      
      ranges <- unlist(ranges(ranges), use.names=FALSE)
      .getExons(txdb=ann, ranges=ranges,
                chrom=chrom, strand=strand,
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
.getExons <- function(txdb, ranges=NULL, chrom=NULL,
                      strand=NULL,
                      rangeRestr=c("both","either","start","end"),
                      expand=FALSE) {
  rangeRestr <- match.arg(rangeRestr)
  len <- max(length(chrom), length(strand), length(ranges))
  ## Note that .mapExons() uses the same SQL query.
  sql <- paste("SELECT exon_id, exon_chrom, exon_strand,",
               "exon_start, exon_end",
               "FROM exon INNER JOIN exon_rtree",
               "ON (exon._exon_id=exon_rtree._exon_id")
  if (len > 0) {
    if (!is.null(chrom)) {
      if (length(chrom) < len)
        chrom <- rep(chrom, length.out = len)
      sqlwhere <- paste("exon_chrom='", chrom, "'", sep="")
    } else {
      sqlwhere <- character(0)
    }
    if (!is.null(strand)) {
      if (length(strand) < len)
        strand <- rep(strand, length.out = len)
      sqladd <- paste("exon_strand='", strand, "'", sep="")
      if (length(sqlwhere) == 0)
        sqlwhere <- sqladd
      else
        sqlwhere <- paste(sqlwhere, sqladd, sep = " AND ")
    }
    if (!is.null(ranges)) {
      if (length(ranges) < len)
        ranges <- rep(ranges, length.out = len)
      sqladd <- .makeRangeSQL("exon_start", "exon_end", ranges, rangeRestr)
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
  sql <- paste(sql, "ORDER BY exon_id")
  if (getOption("verbose", FALSE)) {
    .printSQL(sql)
  }
  ans <- dbGetQuery(txdb@conn, sql)
  ans <- RangedData(ranges = IRanges(start = ans[["exon_start"]],
                                     end   = ans[["exon_end"]]),
                    strand = strand(ans[["exon_strand"]]),
                    space  = ans[["exon_chrom"]],
                    GF_exonId = ans[["exon_id"]])
  if (expand) {
    sqltx <- paste("SELECT exon_id, tx_id",
                   "FROM exon",
                   "INNER JOIN splicing",
                   "ON (exon._exon_id=splicing._exon_id)",
                   "INNER JOIN transcript",
                   "ON (splicing._tx_id=transcript._tx_id)")
    if (len != 0)
      sqltx <- paste(sqltx,
                     "INNER JOIN exon_rtree",
                     "ON (exon._exon_id=exon_rtree._exon_id)",
                     "WHERE", sqlwhere)
    sqltx <- paste(sqltx, "ORDER BY exon_id")
    if (getOption("verbose", FALSE)) {
      .printSQL(sqltx)
    }
    tx <- dbGetQuery(txdb@conn, sqltx)
    ans[["transcript"]] <-
      IRanges:::newCompressedList("CompressedCharacterList",
                                  as.character(tx[["tx_id"]]),
                                  end = end(Rle(tx[["exon_id"]])))
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

