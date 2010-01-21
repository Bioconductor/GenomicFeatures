

transcriptsByRanges <- function(txdb, ranges, restrict = "any",
                                columns=c("tx_id", "tx_name"))
{
  ## check/assign the strand, ranges c-somes
  if(.checkFields(txdb, "tx_strand", ranges[["strand"]], "transcript")){
    strand <- ranges[["strand"]]
  }else{stop("Strand values for ranges do not match annotation DB")}
  if(.checkFields(txdb, "tx_chrom", space(ranges), "transcript")){
    chrom <- space(ranges)
  }else{stop("Space values for ranges do not match annotation DB")}
  ranges <- unlist(ranges(ranges), use.names=FALSE)
  .mapTranscripts(txdb=txdb, ranges=ranges,
                  chrom=chrom, strand=strand,
                  type = restrict, col=columns,
                  format="get")
}



## ## If there is not ranged Data object, then we just want it all...
## setMethod("transcriptsByRanges", "missing",
##     function(txdb, ranges, restrict = "any", columns=c("tx_id", "tx_name"))
##     {
##       .mapTranscripts(txdb=txdb, ranges=NULL,
##                       chrom=NULL, strand=NULL,
##                       type=restrict, col=columns,
##                       format="get")
##     }
## )


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



## TODO: Need to handle the case where we don't have a ranges... (missing
## method is gone)
exonsByRanges <- function(txdb, ranges, restrict = "any",
                          columns=c("exon_id", "exon_name"))
{
  ## check/assign the strand, ranges c-somes
  if(.checkFields(txdb, "exon_strand", ranges[["exon_strand"]], "exon")){
    strand <- ranges[["exon_strand"]]
  }else{stop("Strand values for ranges do not match annotation DB")}
  if(.checkFields(txdb, "exon_chrom", space(ranges), "exon")){
    chrom <- space(ranges)
  }else{stop("Space values for ranges do not match annotation DB")}
  ranges <- unlist(ranges(ranges), use.names=FALSE)
  .mapExons(txdb=txdb, ranges=ranges,
            chrom=chrom, strand=strand,
            type = restrict, col=columns,
            format="get")
}



## ## If there is not ranged Data object, then we want it all...
## setMethod("exonsByRanges", "missing",
##     function(txdb, ranges, restrict = "any", columns=c("exon_id", "exon_name"))
##     {
##       .mapExons(txdb=txdb, ranges=NULL,
##                 chrom=NULL, strand=NULL,
##                 type=restrict, col=columns,
##                 format="get")
##     }
## )




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






##Time for methods to access data via more direct queries.

.appendValsSQL <- function(vals){
  vals = unlist2(vals) ##AnnotationDbi:::unlist2 preserves names with repeats
  sql <- character()
  for(i in seq_len(length(vals))){
    sql <- c(sql, paste("AND ",names(vals[i])," = ","'",vals[i],"'",sep=""))
  }
  paste(sql, collapse=" ")
}

## This is the core function for looking up transcripts
.lookupTranscripts <- function(txdb, vals, col=c("tx_id", "tx_name")){
  ## check that txdb is in fact a TranscriptDb object
  if(!is(txdb,"TranscriptDb"))stop("txdb MUST be a TranscriptDb object.")
  
  ## check the vals:
  valNames <- c("gene_id", "tx_id", "tx_name", "tx_chrom", "tx_strand")
  if(!all(names(vals) %in% valNames)){
    stop(paste("Argument names for vals must be some combination of: ",
               valNames,sep=""))
  }

  ## check the cols:
  colNames <- c("tx_id", "tx_name", "gene_id", "exon_id","cds_id",
                NULL, character(0))
  if(!all(col %in% colNames)){
    stop(paste("Arguments to column must be some combination of: ",
               colNames,sep=""))
  }

  ## base SQL query 
  sql <- paste("SELECT gene_id, t._tx_id AS tx_id, tx_name, tx_chrom, tx_strand,",
               "tx_start, tx_end, _exon_id AS exon_id, _cds_id AS cds_id",
               "FROM transcript AS t, transcript_rtree AS trt,",
               "gene AS g, splicing AS s",
               "WHERE t._tx_id=trt._tx_id AND t._tx_id=g._tx_id",
               "AND t._tx_id=s._tx_id")


  ## Now we just need to finish the where clause with stuff in "vals"
  sql <- paste(sql, .appendValsSQL(vals), collapse=" ")
  ans <- dbGetQuery(txdb@conn, sql)
  
  if(dim(ans)[1] >0){
      rd <- .formatRD(ans, "get", "tx")
      if(is.null(col) || length(col)==0){
        return(rd)
      }else{
        return(.appendCols(rd, ans, col))
      }
  }else{warning("Please be advised that no matching data was found.")}
}


## 
transcripts <- function(txdb, vals, columns=c("tx_id", "tx_name"))
{
  ##Add error here
  if(is.data.frame(txdb) && is.integer(vals)){stop("This is not the transcripts function that you are looking for. Please use transcripts_deprecated instead.")}
  .lookupTranscripts(txdb=txdb, vals, col=columns)
}
