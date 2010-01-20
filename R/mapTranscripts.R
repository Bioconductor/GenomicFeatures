##This should be similar to getTranscripts, but is for the more commonly used
##mapTranscripts() and mapExons() methods.



## For looping we want to get data for just the unique chrom/strand
## combinations that actually occur...
## BUT if one or both of the chrom/strand is NULL: we must expand
.createUniqueChromStrand <- function(txdb, table, chrom, strand, colnames){
  if(is.null(chrom) && !is.null(strand)){
    allChrms = .getValsFromTable(txdb, colnames[1L], table)
    uChromStrand <- unique(cbind(rep(allChrms,length(strand)),
                                 rep(strand,each=length(allChrms))) )  
  }
  if(!is.null(chrom) && is.null(strand)){
    allStrnds = .getValsFromTable(txdb, colnames[2L], table)
    uChromStrand <- unique(cbind(rep(chrom,length(allStrnds)),
                                 rep(allStrnds,each=length(chrom))) )  
  }
  if(is.null(strand) && is.null(chrom)){
    allChrms = .getValsFromTable(txdb, colnames[1L], table)
    allStrnds = .getValsFromTable(txdb, colnames[2L], table)
    uChromStrand <- unique(cbind(rep(allChrms,length(allStrnds)),
                                 rep(allStrnds,each=length(allChrms))))
  }
  if(!is.null(chrom) && !is.null(strand)){
#    uChromStrand <- unique(cbind(chrom,strand))
    uChromStrand <- unique(cbind(as.character(chrom),as.character(strand)))
  }
  uChromStrand
}


.setRange <- function(ranges, txdb, chrom, strand, chr, strnd){
  if(is.null(chrom) && !is.null(strand)){
    range <- ranges[strand==strnd]
  }
  if(!is.null(chrom) && is.null(strand)){
    range <- ranges[chrom==chr]
  }
  if(is.null(strand) && is.null(chrom)){
    range <- ranges
  }
  if(!is.null(chrom) && !is.null(strand)){
    range <- ranges[chrom==chr & strand==strnd]
  }
  range
}


.mapFeatures <- function(txdb, sql, range, chrom, strnd, colprefix,
                         type="any"){

  sqlChrom <- paste(" AND ",colprefix,"_chrom = '",chrom,"'",sep="")
  sqlStrand <- paste(" AND ",colprefix,"_strand = '",strnd,"'",sep="")
  sql2 <- paste(sql, sqlChrom, sqlStrand)
  res <- dbGetQuery(txdb@conn, sql2)
  
  ## Ranges is the query, the DB returned tx_start and tx_end will be the
  ## subject
  startName <- paste(colprefix,"_start",sep="")
  endName <- paste(colprefix,"_end",sep="")
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


.appendCols <- function(rd, ans, col){## for each value in cols
  ## verify that cols are in ans
  if(all(col %in% colnames(ans))){
    for(i in seq_len(length(col))){
      rd[[col[i]]] <- ans[[col[i]]]
    }
  }else{stop(paste("There is either a problem with the requested columns",
                   "or the answer returned from the DB is not complete."))}
  rd
}



.mapTranscripts <- function(txdb, ranges=NULL, chrom=NULL,
                            strand=NULL, type="any",
                            col=c("tx_id", "tx_name")) {
  ## check the cols:
  colNames <- c("tx_id", "tx_name", "gene_id", "exon_id","cds_id",
                NULL, character(0))
  if(!all(col %in% colNames)){
    stop(paste("Arguments to column must be some combination of: ",
               colNames,sep=""))
  }
  ## Note that .getTranscripts() uses the same SQL query.
  sql <- paste("SELECT gene_id, t._tx_id AS tx_id, tx_name, tx_chrom, tx_strand,",
               "tx_start, tx_end, _exon_id AS exon_id, _cds_id AS cds_id",
               "FROM transcript AS t, transcript_rtree AS trt,",
               "gene AS g, splicing AS s",
               "WHERE t._tx_id=trt._tx_id AND t._tx_id=g._tx_id",
               "AND t._tx_id=s._tx_id")
  ## We could speed this up by doing less complicated joins when possible...
  ## The length of the range, the chromosome and the strand must be the same.
  len <- max(length(chrom), length(strand), length(ranges))
  if(len <= 0){stop("Ranges, chromosomes and strands are all required.")}
  if (!is.null(chrom) && length(chrom) < len){
    stop("Not enough chromosomes for all the ranges")}
  if (!is.null(strand) && length(strand) < len){
    stop("Not enough strands for all the ranges")}
  
  uChromStrand <- .createUniqueChromStrand(txdb,"transcript",
                                           chrom, strand,
                                           c("tx_chrom", "tx_strand"))
  ans <- data.frame()
  for(i in seq_len(dim(uChromStrand)[1])){
    chr <- uChromStrand[i,1]
    strnd <- uChromStrand[i,2]
    range <- .setRange(ranges, txdb, chrom, strand, chr, strnd)
    if(length(range) > 0){
      newData <- .mapFeatures(txdb, sql, range, chr, strnd, "tx", type)
      if(dim(newData)[1] > 0){
        chrs <- rep(chr, dim(newData)[1])
        strnds <-rep(strnd, dim(newData)[1])
        ans <- rbind(ans, cbind(chrs,strnds,newData))
      }
    }
  }
  if(dim(ans)[1] >0){
      rd =  RangedData(ranges = IRanges(start = ans[["start"]],
                          end = ans[["end"]]),
                       strand = ans[["tx_strand"]],
                        space = ans[["tx_chrom"]])
      if(is.null(col) || length(col)==0){
        return(rd)
      }else{
        return(.appendCols(rd, ans, col))
      }
  }else{warning("Please be advised that no matching data was found.")}
}




setMethod("transcriptsByRanges", "RangedData",
    function(txdb, ranges, restrict = "any", columns=c("tx_id", "tx_name"))
    {
      ## check that txdb is in fact a TranscriptDb object
      if(!is(txdb,"TranscriptDb"))stop("txdb MUST be a TranscriptDb object.")
      ## check that any supplied strands are what we expect.
      ## note that NULL is ok and will not trip this error.
      ## If there are NULL strands or chromosomes, then .mapTranscripts() will
      ## have to search a bit differently
      if(.checkFields(txdb, "tx_strand", ranges[["strand"]], "transcript")){
        strand <- ranges[["strand"]]
      }else{stop("Strand values for ranges do not match annotation DB")}
      ## check that the chromosomes are what we expect. 
      if(.checkFields(txdb, "tx_chrom", space(ranges), "transcript")){
        chrom <- space(ranges)
      }else{stop("Space values for ranges do not match annotation DB")}
      ranges <- unlist(ranges(ranges), use.names=FALSE)
      .mapTranscripts(txdb=txdb, ranges=ranges,
                      chrom=chrom, strand=strand,
                      type = restrict, col=columns)
    }
)



## ## Test code:
## library(GenomicFeatures)
## library(DBI)
## library(IRanges)
## ranges <- IRanges(start = c(1,20,300,254872),end =c(22,41,321,254872+21))
## chrom <- c("chr1", "chr1", "chr2", "chr2")
## strand <- c("+", "-", "+", "+")
## txdb <- loadFeatures("testDB.sqlite")
## #GenomicFeatures:::.mapTranscripts(txdb, ranges, chrom, strand)

## #rd = RangedData(ranges, space=chrom, strand=strand)
## #mapTranscripts(rd, txdb)










##Additional code to support mapping on exons:


.mapExons <- function(txdb, ranges=NULL, chrom=NULL,
                      strand=NULL, type="any",
                      col=c("exon_id")) {
  ## check the cols:
  colNames=c("exon_id", NULL, character(0))
  if(!all(col %in% colNames)){
    stop(paste("Arguments to column must be some combination of: ",
               colNames,sep=""))
  }  
  ## Note that .getExons() uses the same SQL query.
  sql <- paste("SELECT exon_id, exon_chrom, exon_strand,",
               "exon_start, exon_end",
               "FROM exon INNER JOIN exon_rtree",
               "ON (exon._exon_id=exon_rtree._exon_id")
  ## The length of the range, the chromosome and the strand must be the same.
  len <- max(length(chrom), length(strand), length(ranges))
  if(len <= 0){stop("Ranges, chrom and strands are all required.")}
  if (!is.null(chrom) && length(chrom) < len){
    stop("Not enough chromosomes for all the ranges")}
  if (!is.null(strand) && length(strand) < len){
    stop("Not enough strands for all the ranges")}
  uChromStrand <- .createUniqueChromStrand(txdb,"exon",
                                           chrom, strand,
                                           c("exon_chrom", "exon_strand"))
  ans <- data.frame()
  for(i in seq_len(dim(uChromStrand)[1])){
    chr <- uChromStrand[i,1]
    strnd <- uChromStrand[i,2]
    range <- .setRange(ranges, txdb, chrom, strand, chr, strnd)
    if(length(range) > 0){
      newData <- .mapFeatures(txdb, sql, range, chr, strnd , "exon", type)
      if(dim(newData)[1] > 0){
        chrs <- rep(chr, dim(newData)[1])
        strnds <-rep(strnd, dim(newData)[1])
        ans <- rbind(ans, cbind(chrs,strnds,newData))
      }
    }
  }
  if(dim(ans)[1] >0){
      rd =  RangedData(ranges = IRanges(start = ans[["start"]],
                          end = ans[["end"]]),
                       strand = ans[["exon_strand"]],
                        space = ans[["exon_chrom"]])
      if(is.null(col) || length(col)==0){
        return(rd)
      }else{
        return(.appendCols(rd, ans, col))
      }    
    }else{warning("Please be advised that no matching data was found.")}
}


setMethod("exonsByRanges", "RangedData",
    function(txdb, ranges, restrict = "any", columns=c("exon_id", "exon_name"))
    {
      ## check that txdb is in fact a TranscriptDb object
      if(!is(txdb,"TranscriptDb"))stop("txdb MUST be a TranscriptDb object.")
      ## check the strand and c-somes
      if(.checkFields(txdb, "exon_strand", ranges[["strand"]], "exon")){
        strand <- ranges[["strand"]]
      }else{stop("Strand values for ranges do not match annotation DB")}
      if(.checkFields(txdb, "exon_chrom", space(ranges), "exon")){
        chrom <- space(ranges)
      }else{stop("Space values for ranges do not match annotation DB")}
      ranges <- unlist(ranges(ranges), use.names=FALSE)
      .mapExons(txdb=txdb, ranges=ranges,
                chrom=chrom, strand=strand,
                type = restrict, col=columns)
    }
)


## Good news: this seems to work generically as expected.

## ## Test code:
## library(GenomicFeatures)
## library(DBI)
## library(IRanges)
## ranges <- IRanges(start = c(1,20,300,254872),end =c(22,41,321,254872+21))
## chrom <- c("chr1", "chr1", "chr2", "chr2")
## strand <- c("+", "-", "+", "+")
## txdb <- loadFeatures("testDB.sqlite")
## #GenomicFeatures:::.mapExons(txdb, ranges, chrom, strand)

## So now I just have to wrap this in a method and make it into a nice example.
## rd = RangedData(ranges, space = chrom, strand = strand)
## mapExons(rd, txdb)


## TODO: edit findOverlaps so that I can use a different type parameter...
## type can be:  type = c("any", "start", "finish", "during", "equal")




## reproduce Patricks bug
## library(GenomicFeatures)
## library(DBI)
## library(IRanges)
## ranges <- IRanges(start = c(1,20,300,254872),end =c(22,41,321,254872+21))
## chrom <- c("chr1", "chr1", "chr2", "chr2")
## strand <- c("+", "-", "+", "+")
## txdb <- loadFeatures("testDB.sqlite")

## system.time(getTranscripts(txdb=txdb))
## rd = getTranscripts(txdb=txdb)
## system.time(mapTranscripts(rd,txdb))
