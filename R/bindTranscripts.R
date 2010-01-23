##TEMP (just untill we refactor bindTranscripts and bindExons)
.mapTranscripts <- function(txdb, ranges=NULL, chrom=NULL,
                            strand=NULL, type="any",
                            col=c("tx_id", "tx_name"),
                            format=c("map","get")) {
  ## check that txdb is in fact a TranscriptDb object
  if(!is(txdb,"TranscriptDb"))stop("txdb MUST be a TranscriptDb object.")
  ## check the format:
  format=match.arg(format)
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
      rd <- .formatRD(ans, format, "tx")
      if(is.null(col) || length(col)==0){
        return(rd)
      }else{
        return(.appendCols(rd, ans, col))
      }
  }else{warning("Please be advised that no matching data was found.")}
}


## This is the core function for mapping transcripts
.mapExons <- function(txdb, ranges=NULL, chrom=NULL,
                      strand=NULL, type="any",
                      col=c("exon_id"),
                      format=c("map","get")) {
  ## check that txdb is in fact a TranscriptDb object
  if(!is(txdb,"TranscriptDb"))stop("txdb MUST be a TranscriptDb object.")
  ## check the format:
  format=match.arg(format)
  ## check the cols:
  colNames=c("exon_id", NULL, character(0))
  if(!all(col %in% colNames)){
    stop(paste("Arguments to column must be some combination of: ",
               colNames,sep=""))
  }  
  ## base SQL query:
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
      rd <- .formatRD(ans, format, "exon")
      if(is.null(col) || length(col)==0){
        return(rd)
      }else{
        return(.appendCols(rd, ans, col))
      }    
    }else{warning("Please be advised that no matching data was found.")}
}

##END TEMP

bindTranscripts <- function(txdb, ranges, restrict = "any",
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
                  type=restrict, col=columns,
                  format="map")
}




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




bindExons <- function(txdb, ranges, restrict = "any",
                      columns=c("exon_id", "exon_name"))
{
  ## check/assign the strand, ranges c-somes
  if(.checkFields(txdb, "exon_strand", ranges[["strand"]], "exon")){
    strand <- ranges[["strand"]]
  }else{stop("Strand values for ranges do not match annotation DB")}
  if(.checkFields(txdb, "exon_chrom", space(ranges), "exon")){
    chrom <- space(ranges)
  }else{stop("Space values for ranges do not match annotation DB")}
  ranges <- unlist(ranges(ranges), use.names=FALSE)
  .mapExons(txdb=txdb, ranges=ranges,
            chrom=chrom, strand=strand,
            type = restrict, col=columns,
            format="map")
}



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


