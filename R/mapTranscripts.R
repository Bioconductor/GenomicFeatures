##This should be similar to getTranscripts, but is for the more commonly used
##mapTranscripts() and mapExons() methods.



## For looping we want to get data for just the unique chromosome and strand
## combinations that actually occur...
## BUT if one or both of the chromosome/strand is NULL: we must expand
.createUniqueChromStrand <- function(ann, table, chromosome, strand){
  if(is.null(chromosome) && !is.null(strand)){
    allChrms = .getValsFromTable(ann, "chromosome", table)
    uChromStrand <- unique(cbind(rep(allChrms,length(strand)),
                                 rep(strand,each=length(allChrms))) )  
  }
  if(!is.null(chromosome) && is.null(strand)){
    allStrnds = .getValsFromTable(ann, "strand", table)
    uChromStrand <- unique(cbind(rep(chromosome,length(allStrnds)),
                                 rep(allStrnds,each=length(chromosome))) )  
  }
  if(is.null(strand) && is.null(chromosome)){
    allChrms = .getValsFromTable(ann, "chromosome", table)
    allStrnds = .getValsFromTable(ann, "strand", table)
    uChromStrand <- unique(cbind(rep(allChrms,length(allStrnds)),
                                 rep(allStrnds,each=length(allChrms))))
  }
  if(!is.null(chromosome) && !is.null(strand)){
#    uChromStrand <- unique(cbind(chromosome,strand))
    uChromStrand <- unique(cbind(as.character(chromosome),as.character(strand)))
  }
  uChromStrand
}


.setRange <- function(ranges, ann, chromosome, strand, chrom, strnd){
  if(is.null(chromosome) && !is.null(strand)){
    range <- ranges[strand==strnd]
  }
  if(!is.null(chromosome) && is.null(strand)){
    range <- ranges[chromosome==chrom]
  }
  if(is.null(strand) && is.null(chromosome)){
    range <- ranges
  }
  if(!is.null(chromosome) && !is.null(strand)){
    range <- ranges[chromosome==chrom & strand==strnd]
  }
  range
}


.mapFeatures <- function(ann, sql, range, chrom, strnd, abr, fldAbr,
                         type="any"){

  sqlChrom <- paste(" AND ",abr,".chromosome = '",chrom,"'",sep="")
  sqlStrand <- paste(" AND ",abr,".strand = '",strnd,"'",sep="")
  sql2 <- paste(sql, sqlChrom, sqlStrand)
  res <- dbGetQuery(ann@conn, sql2)
  
  ## Ranges is the query, the DB returned tx_start and tx_end will be the
  ## subject
  startName <- paste(fldAbr,"_start",sep="")
  endName <- paste(fldAbr,"_end",sep="")
  if(dim(res)[1] > 0){
    map <- findOverlaps(range,
                        IRanges(start = unlist(res[startName]),
                                end = unlist(res[endName])),
                        multiple=TRUE,
                        type=type)
  
    ##Then subset the query result with only the matching stuff.
    ##And I also have to match that up with the stuff in range...
    matchRange <- range[map@matchMatrix[,1]]
    if(!is.na(map@matchMatrix[1])){
      ans <- cbind(as.data.frame(matchRange),res[map@matchMatrix[,2],])
    }else ans <- data.frame()
  }
  ans
}


.mapTranscripts <- function(ann, ranges=NULL, chromosome=NULL,
                            strand=NULL, type="any") {
  len <- max(length(chromosome), length(strand), length(ranges))
  sql <- paste("SELECT t.tx_id, t.tx_name, t.chromosome, t.strand,",
               "tt.tx_start, tt.tx_end",
               "FROM transcripts AS t,",
               "transcripts_rtree AS tt",
               "WHERE t._tx_id=tt._tx_id")  
  ## The length of the range, the chromosome and the strand must be the same.
  len <- max(length(chromosome), length(strand), length(ranges))
  if(len <= 0){stop("Ranges, chromosomes and strands are all required.")}
  if (!is.null(chromosome) && length(chromosome) < len){
    stop("Not enough chromosomes for all the ranges")}
  if (!is.null(strand) && length(strand) < len){
    stop("Not enough strands for all the ranges")}
  
  uChromStrand <- .createUniqueChromStrand(ann,"transcripts",
                                           chromosome, strand)
  ans <- data.frame()
  for(i in seq_len(dim(uChromStrand)[1])){
    chrom <- uChromStrand[i,1]
    strnd <- uChromStrand[i,2]
    range <- .setRange(ranges, ann, chromosome, strand, chrom, strnd)
    if(length(range) > 0){
      newData <- .mapFeatures(ann, sql, range, chrom, strnd ,"t", "tx", type)
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
               GF_txId    = ans[["tx_id"]],
##                txStart = ans[["tx_end"]],
##                txEnd = ans[["tx_start"]],
               txName = ans[["tx_name"]],##  temp.
               space      = ans[["chromosome"]])
  }else{warning("Please be advised that no matching data was found.")}
}



setMethod("mapTranscripts", "RangedData",
    function(ranges, ann, rangeRestr = "any")
    {
      ## check that any supplied strands are what we expect.
      ## note that NULL is ok and will not trip this error.
      ## If there are NULL strands or chromosomes, then .mapTranscripts() will
      ## have to search a bit differently
      if(.checkFields(ann, "strand", ranges[["strand"]], "transcripts")){
        strand <- ranges[["strand"]]
      }else{stop("Strand values for ranges do not match annotation DB")}
      ## check that the chromosomes are what we expect. 
      if(.checkFields(ann, "chromosome", space(ranges), "transcripts")){
        chromosome <- space(ranges)
      }else{stop("Space values for ranges do not match annotation DB")}
      ranges <- unlist(ranges(ranges), use.names=FALSE)
      .mapTranscripts(ann=ann, ranges=ranges,
                      chromosome=chromosome, strand=strand,
                      type = rangeRestr)
    }
)



## ## Test code:
## library(GenomicFeatures)
## library(DBI)
## library(IRanges)
## ranges <- IRanges(start = c(1,20,300,254872),end =c(22,41,321,254872+21))
## chromosome <- c("chr1", "chr1", "chr2", "chr2")
## strand <- c("+", "-", "+", "+")
## ann <- loadFeatures("testDB.sqlite")
## #GenomicFeatures:::.mapTranscripts(ann, ranges, chromosome, strand)

## #rd = RangedData(ranges, space=chromosome, strand=strand)
## #mapTranscripts(rd, ann)






##Additional code to support mapping on exons:


.mapExons <- function(ann, ranges=NULL, chromosome=NULL,
                            strand=NULL, type="any") {
  len <- max(length(chromosome), length(strand), length(ranges))
  sql <- paste("SELECT e.exon_id, e.chromosome, e.strand,",
               "ee.exon_start, ee.exon_end",
               "FROM exons AS e,",
               "exons_rtree AS ee",
               "WHERE e._exon_id=ee._exon_id")
  ## The length of the range, the chromosome and the strand must be the same.
  len <- max(length(chromosome), length(strand), length(ranges))
  if(len <= 0){stop("Ranges, chromosomes and strands are all required.")}
  if (!is.null(chromosome) && length(chromosome) < len){
    stop("Not enough chromosomes for all the ranges")}
  if (!is.null(strand) && length(strand) < len){
    stop("Not enough strands for all the ranges")}
  uChromStrand <- .createUniqueChromStrand(ann,"exons",
                                           chromosome, strand)
  ans <- data.frame()
  for(i in seq_len(dim(uChromStrand)[1])){
    chrom <- uChromStrand[i,1]
    strnd <- uChromStrand[i,2]
    range <- .setRange(ranges, ann, chromosome, strand, chrom, strnd)
    if(length(range) > 0){
      newData <- .mapFeatures(ann, sql, range, chrom, strnd ,"e", "exon", type)
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
               GF_exonId  = ans[["exon_id"]],
##                exonStart = ans[["exon_start"]],
##                exonEnd = ans[["exon_end"]],
               space      = ans[["chromosome"]])
  }else{warning("Please be advised that no matching data was found.")}
}


setMethod("mapExons", "RangedData",
    function(ranges, ann, rangeRestr = "any")
    {
      ## check the strand and c-somes
      if(.checkFields(ann, "strand", ranges[["strand"]], "exons")){
        strand <- ranges[["strand"]]
      }else{stop("Strand values for ranges do not match annotation DB")}
      if(.checkFields(ann, "chromosome", space(ranges), "exons")){
        chromosome <- space(ranges)
      }else{stop("Space values for ranges do not match annotation DB")}
      ranges <- unlist(ranges(ranges), use.names=FALSE)
      .mapExons(ann=ann, ranges=ranges,
                chromosome=chromosome, strand=strand,
                type = rangeRestr)
    }
)


## Good news: this seems to work generically as expected.

## ## Test code:
## library(GenomicFeatures)
## library(DBI)
## library(IRanges)
## ranges <- IRanges(start = c(1,20,300,254872),end =c(22,41,321,254872+21))
## chromosome <- c("chr1", "chr1", "chr2", "chr2")
## strand <- c("+", "-", "+", "+")
## ann <- loadFeatures("testDB.sqlite")
## #GenomicFeatures:::.mapExons(ann, ranges, chromosome, strand)

## So now I just have to wrap this in a method and make it into a nice example.
## rd = RangedData(ranges, space = chromosome, strand = strand)
## mapExons(rd, ann)


## TODO: edit findOverlaps so that I can use a different type parameter...
## type can be:  type = c("any", "start", "finish", "during", "equal")




## reproduce Patricks bug
## library(GenomicFeatures)
## library(DBI)
## library(IRanges)
## ranges <- IRanges(start = c(1,20,300,254872),end =c(22,41,321,254872+21))
## chromosome <- c("chr1", "chr1", "chr2", "chr2")
## strand <- c("+", "-", "+", "+")
## ann <- loadFeatures("testDB.sqlite")

## system.time(getTranscripts(ann=ann))
## rd = getTranscripts(ann=ann)
## system.time(mapTranscripts(rd,ann))
