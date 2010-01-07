##This should be similar to getTranscripts, but is for the more commonly used
##mapTranscripts() and mapExons() methods.


## Get the list of possible values for a type from a table:
.getValsFromTable <- function(txann, field, table){
  sql <- paste("SELECT",field,"FROM",table)
  as.character(unlist(unique(dbGetQuery(txann@conn, sql))))
}

## Check that the type of thing (chromosome, strand etc.) is in the
## transcripts table
.checkTxFields <- function(txann, type = c("chromosome","strand"), data){
  type <- match.arg(type)
  annot <- .getValsFromTable(txann, type, "transcripts") 
  all(data %in% annot)
}



## For looping we want to get data for just the unique chromosome and strand
## combinations that actually occur...
## BUT if one or both of the chromosome/strand is NULL: we must expand
.createUniqueChromStrand <- function(txann, table, chromosome, strand){
  if(is.null(chromosome) && !is.null(strand)){
    allChrms = .getValsFromTable(txann, "chromosome", table)
    uChromStrand <- unique(cbind(rep(allChrms,length(strand)),
                                 rep(strand,each=length(allChrms))) )  
  }
  if(!is.null(chromosome) && is.null(strand)){
    allStrnds = .getValsFromTable(txann, "strand", table)
    uChromStrand <- unique(cbind(rep(chromosome,length(allStrnds)),
                                 rep(allStrnds,each=length(chromosome))) )  
  }
  if(is.null(strand) && is.null(chromosome)){
    allChrms = .getValsFromTable(txann, "chromosome", table)
    allStrnds = .getValsFromTable(txann, "strand", table)
    uChromStrand <- unique(cbind(rep(allChrms,length(allStrnds)),
                                 rep(allStrnds,each=length(allChrms))))
  }
  if(!is.null(chromosome) && !is.null(strand)){
    uChromStrand <- unique(cbind(chromosome,strand))
  }
  uChromStrand
}



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
  sql <- paste("SELECT t._tx_id, t.chromosome, t.strand,",
               "tt.tx_start, tt.tx_end",
               "FROM transcripts AS t,",
               "transcripts_rtree AS tt",
               "WHERE t._tx_id=tt._tx_id")
## I can't use a rangeRestr parameter here because we are using findOverlaps
## instead of the DB. :(
  
## The length of the range, the chromosome and the strand must be the same.
  len <- max(length(chromosome), length(strand), length(ranges))
  if(len <= 0){stop("Ranges, chromosomes and strands are all required.")}
  if (!is.null(chromosome) && length(chromosome) < len){
    stop("Not enough chromosomes for all the ranges")}
  if (!is.null(strand) && length(strand) < len){
    stop("Not enough strands for all the ranges")}
  
  uChromStrand <- .createUniqueChromStrand(txann,"transcripts",
                                           chromosome, strand)

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

  if(dim(ans)[1] >0){
    RangedData(ranges     = IRanges(start = ans[["start"]],
                                    end   = ans[["end"]]),
               strand     = ans[["strand"]],
               GF_txId    = ans[["_tx_id"]],
               space      = ans[["chromosome"]])
  }else{stop("No matching data found.")}
  
}



setMethod("mapTranscripts", "RangedData",
    function(ranges, ann)
    {
      ## check that any supplied strands are what we expect.
      ## note that NULL is ok and will not trip this error.
      ## If there are NULL strands or chromosomes, then .mapTranscripts() will
      ## have to search a bit differently
      if(.checkTxFields(ann, "strand", ranges[["strand"]])){
        strand <- ranges[["strand"]]
      }else{stop("Strand values for ranges do not match annotation DB")}

      ## check that the chromosomes are what we expect. 
      if(.checkTxFields(ann, "chromosome", space(ranges))){
        chromosome <- space(ranges)
      }else{stop("Space values for ranges do not match annotation DB")}

      ranges <- unlist(ranges(ranges), use.names=FALSE)
      .mapTranscripts(txann=ann, ranges=ranges,
                      chromosome=chromosome, strand=strand)
    }
)



## ## Test code:
## library(GenomicFeatures)
## library(DBI)
## library(IRanges)
## ranges <- IRanges(start = c(1,20,300,254872),end =c(22,41,321,254872+21))
## chromosome <- c("chr1", "chr1", "chr2", "chr2")
## strand <- c("+", "-", "+", "+")
## txAnn <- loadFeatures("testDB.sqlite")
## #GenomicFeatures:::.mapTranscripts(txAnn, ranges, chromosome, strand)

## #rd = RangedData(ranges, space=chromosome, strand=strand)
## #mapTranscripts(rd, txAnn)

## #rd = RangedData(ranges, space=chromosome)
## #mapTranscripts(rd, txAnn)




## NEXT:
## 1) wrap this in a single method to support:
## input = RangedData, TxAnn ; output = RangedData (with one more col).
## 2) write test code (for checking helper methods)
## 3) write a manual page
## 4) since we decided to only return the tx_id, we might want to change to returning the _tx_id starting right now. - done?
