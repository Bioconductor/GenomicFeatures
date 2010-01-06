##This should be similar to getTranscripts, but is for the more commonly used
##mapTranscripts() and mapExons() methods.

.mapFeatures <- function(txann, sql, range, chrom, strnd){

  sqlChrom <- paste(" AND t.chromosome = '",chrom,"'",sep="")
  sqlStrand <- paste(" AND t.strand = '",strnd,"'",sep="")
  sql2 <- paste(sql, sqlChrom, sqlStrand)
  res <- dbGetQuery(txann@conn, sql2)
  
  ## Ranges is the query, the DB returned tx_start and tx_end will be the
  ## subject
  if(dim(res)[1] > 0){
    map <- findOverlaps(range,
                        IRanges(start = unlist(res["tx_start"]),
                                end = unlist(res["tx_end"])),
                        multiple=TRUE)
  
  
    ##Then subset the query result with only the matching stuff.
    ##And I also have to match that up with the stuff in range...
    matchRange <- range[unique(map@matchMatrix[,1])]
    if(!is.na(map@matchMatrix[1])){
      dupRange <- rep(matchRange, table(map@matchMatrix[,1])[[1]])
      ans <- cbind(as.data.frame(dupRange),res[map@matchMatrix[,2],])
    }else ans <- data.frame()
  }
  ans
}


.mapTranscripts <- function(txann, ranges=NULL, chromosome=NULL,
                            strand=NULL) {
  len <- max(length(chromosome), length(strand), length(ranges))
  sql <- paste("SELECT t.tx_id, t.chromosome, t.strand,",
               "tt.tx_start, tt.tx_end",
               "FROM transcripts AS t,",
               "transcripts_rtree AS tt",
               "WHERE t._tx_id=tt._tx_id")
## I can't use a rangeRestr parameter here because we are using findOverlaps
## instead of the DB. :(
  
## The length of the range, the chromosome and the strand must be the same.
  len <- max(length(chromosome), length(strand), length(ranges))
  if(len <= 0){stop("Ranges, chromosomes and strands are all required.")}
  if (length(chromosome) < len){
    chromosome <- rep(chromosome, length.out = len)}
  if (length(strand) < len){
    strand <- rep(strand, length.out = len)}

## For looping we want to get data for just the unique chromosome and strand
## combinations that actually occur...
  uChromStrand <- unique(cbind(chromosome,strand))
  ans <- data.frame()
  for(i in seq_len(dim(uChromStrand)[1])){
    chrom <- uChromStrand[i,1]
    strnd <- uChromStrand[i,2]
    range <- ranges[chromosome==chrom & strand==strnd]
    if(length(range) > 0){
      newData <- .mapFeatures(txann, sql, range, chrom, strnd)
      if(dim(newData)[1] > 0){
        chrs <- rep(chrom, dim(newData)[1])
        strnds <-rep(strnd, dim(newData)[1])
        ans <- rbind(ans, cbind(chrs,strnds,newData))
      }
    }
  }
  
  RangedData(ranges     = IRanges(start = ans[["start"]],
                                  end   = ans[["end"]]),
             strand     = strand(ans[["strand"]]),
             transcript = ans[["tx_id"]],
             space      = ans[["chromosome"]],
             tx_start   = ans[["tx_start"]],
             tx_start   = ans[["tx_end"]])
}






##Test code:
## library(GenomicFeatures)
## library(DBI)
## library(IRanges)
## ranges <- IRanges(start = c(1,20,300,254872),end =c(22,41,321,254872+21))
## chromosome <- c("chr1", "chr1", "chr2", "chr2")
## strand <- c("+", "-", "+", "+")
## txAnn <- loadFeatures("testDB.sqlite")
## GenomicFeatures:::.mapTranscripts(txAnn, ranges, chromosome, strand)




