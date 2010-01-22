
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


