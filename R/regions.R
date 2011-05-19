## 'genes' should be data.frame from dump of UCSC's "knownGenes" table
## Eventually, that should be available from an annotation package
## 'proximal' is the number of bases on either side of TSS and 3'-end for
##   the promoter and end region, respectively.
## 'distal' is the number of bases on either side for upstream/downstream,
##   i.e. enhancer/silencer regions.
transcripts_deprecated <- function(genes, proximal = 500L, distal = 10000L) {
  msg <- c("  transcripts_deprecated() is deprecated.\n",
           "  Please use transcripts() on a TranscriptDb object instead ",
           "(see ?transcripts).")
  .Deprecated(msg=msg)
  transcript <- IRanges(genes$txStart, genes$txEnd)
  # some have multiple exons, CDS, etc
  notdup <- !duplicated(genes[c("chrom", "txStart", "txEnd")])
  transcript <- transcript[notdup]

  ## direction of transcription depends on strand (but UCSC is only by +)
  pos <- genes$strand[notdup] == "+"
  starts <- ifelse(pos, start(transcript), end(transcript))
  ends <- ifelse(pos, end(transcript), start(transcript))
  proximal <- ifelse(pos, proximal, -proximal)
  distal <- ifelse(pos, distal, -distal)
  offset <- ifelse(pos, 1L, -1L)

  rangebind <- function(x, y) IRanges(pmin.int(x,y), pmax.int(x,y))

  ## [start-proximal,start+proximal-1]
  promoter <- rangebind(starts - proximal, starts + proximal - offset)
  ## [end-proximal+1,end+proximal]
  threeprime <- rangebind(ends - proximal + offset, ends + proximal)
  ## [start-distal,start-proximal-1]
  upstream <- rangebind(starts - distal, starts - proximal - offset)
  ## [end+proximal+1, end+distal]
  downstream <- rangebind(ends + proximal + offset, ends + distal)

  RangedData(transcript,
             gene = genes$name[notdup],
             promoter = promoter,
             threeprime = threeprime,
             upstream = upstream,
             downstream = downstream,
             space = genes$chrom[notdup])
}


.countPos <- function(pos) {
   strsplit(as.character(pos), ",", fixed=TRUE)
}

.splitPos <- function(pos) {
    as.integer(unlist(.countPos(pos), use.names=FALSE))
}


## Introns are stored in the local DB ... for now

exons_deprecated <- function(genes) {
  msg <- c("  exons_deprecated() is deprecated.\n",
           "  Please use exons() on a TranscriptDb object instead ",
           "(see ?exons).")
  .Deprecated(msg=msg)
  ## [exon:start, exon:end]
  RangedData(IRanges(.splitPos(genes$exonStarts),
                     .splitPos(genes$exonEnds)),
             gene = rep(genes$name, genes$exonCount),
             space = rep(genes$chrom, genes$exonCount))
}


## RG thinks that there is always one more exon than intron
## and that the first intron comes after the first exon
## For this data at least, I believe that RG is correct. -MC
introns_deprecated <- function(genes) {
  msg <- c("  introns_deprecated() is deprecated.\n",
           "  Please use intronsByTranscript() on a TranscriptDb object instead ",
           "(see ?intronsByTranscript).")
  .Deprecated(msg=msg)
  ## [intron:start, intron:end]
  tmp <- .countPos(genes$intronStarts)
  nIntrons <- sapply(tmp, length)
  RangedData(IRanges(.splitPos(genes$intronStarts),
                     .splitPos(genes$intronEnds)),
             gene = rep(genes$name, nIntrons),
             space = rep(genes$chrom,nIntrons))  
}
