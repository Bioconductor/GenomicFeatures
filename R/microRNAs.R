### =========================================================================
### Extractors for features in other databases
### -------------------------------------------------------------------------
###
### This is for extractors that do NOT point to the TxDb proper.
### Such extractors can point to other databases (mirbase.db) OR they can
### point to other FeatureDbs within the same package.

## helpers for microRNAs

## helpers for for translating chromosomes. Now we have to assume a universal
## translator.  It seems that the chroms are in biomaRt style for mirbase.  So
## for biomaRt, return them as is, but for UCSC, add "chr" prefix.
.translateChromsForUCSC <- function(csomes){
  ## paste0("chr", csomes)
    csomes
}

.translateChromsForBiomaRt <- function(csomes){
  csomes
}

.syncSeqlevel <- function(txdb, ans){
  isActSeq <- .isActiveSeq(txdb)
  n2oNames <- levels(seqnames(ans))
  n2o <- match(seqnames(seqinfo(txdb)), n2oNames)
  seqinfo(ans, new2old=n2o) <- seqinfo(txdb)
  seqlevels(ans, pruning.mode="coarse") <- names(isActSeq)[isActSeq]
  ans
}

## main function
.microRNAs <- function(x){
  ## get the data about whether or not we have any info.
  con <- dbconn(x)
  bld <- dbGetQuery(con,
           "SELECT value FROM metadata WHERE name='miRBase build ID'")
  src <- dbGetQuery(con,
           "SELECT value FROM metadata WHERE name='Data source'")[[1]]

  ## And if not - bail out with message
  if(is.na(bld) || dim(bld)[1]==0){
    stop("this TxDb does not have a miRBase build ID specified")}
  ## now connect to mirbase
  loadNamespace("mirbase.db") ## strictly required

  ## What I need is the join of mirna with mirna_chromosome_build (via _id),
  ## that is then filtered to only have rows that match the species which goes
  ## with the build.
  
  ## connection
  mcon <- mirbase.db::mirbase_dbconn()
  ## 1st lets get the organism abbreviation
  sql <- paste0("SELECT organism FROM mirna_species WHERE genome_assembly",
                " LIKE '", bld, "%'")

  organism <- dbGetQuery(mcon, sql)[[1]]
  ## now get data and make a GRanges from it
  sql <- paste0("SELECT * from mirna_chromosome_build AS csome INNER JOIN ",
                "(SELECT _id,mirna_id,organism from mirna) AS mirna ",
                "WHERE mirna._id=csome._id AND organism='", organism, "' ")
  data <- dbGetQuery(mcon, sql)

  ## convert chromosomes
  csomes <- switch(src,
                   BioMart=.translateChromsForBiomaRt(data$xsome),
                   UCSC=.translateChromsForUCSC(data$xsome),
                   data$xsome)
  ## build GRanges
  ans <- GRanges(seqnames=csomes,
                 ranges=IRanges( ## sign may be reversed
                   start=abs(data$contig_start),
                   end=abs(data$contig_end)),
                 mirna_id = data$mirna_id,
                 strand=data$strand)
  
  ## Filter seqinfo
  .syncSeqlevel(x, ans)
}

setGeneric("microRNAs", function(x) standardGeneric("microRNAs"))

## Then set our method
setMethod("microRNAs", "TxDb", .microRNAs)



## main function
.tRNAs <- function(x) {
  fdbpkg <- "FDb.UCSC.tRNAs"
  fdbenv <- loadNamespace(fdbpkg)
  ## get the current package name
  pkgName <- makePackageName(x)
  ## from here we know what the FDB should MUST look like
  fdbName <- sub("TxDb","FDb",pkgName)
  fdbName <- unlist(strsplit(fdbName,"\\."))
  fdbName[5] <- "tRNAs"
  fdbString <- paste(fdbName,collapse=".")
  if (!exists(fdbString, envir=fdbenv)) {
      stop("there is no tRNA data available for this organism/source")
  } else {
      fdb <- get(fdbString, fdbenv)
      ans <- features(fdb)
  }
  ## Now check active seqs and set the seqlevels
  .syncSeqlevel(x, ans)
}

setGeneric("tRNAs", function(x) standardGeneric("tRNAs"))

setMethod("tRNAs", "TxDb", .tRNAs)


## Test code for new TXTYPE support (BC vs new code)
## library(TxDb.Hsapiens.BioMart.ensembl.GRCh38);txdb2= TxDb.Hsapiens.BioMart.ensembl.GRCh38;transcripts(txdb2, columns='TXTYPE')
## exons(txdb2, columns='TXTYPE')
## And this one works now
## library(TxDb.Hsapiens.UCSC.hg19.knownGene);txdb  = TxDb.Hsapiens.UCSC.hg19.knownGene;transcripts(txdb, columns='TXTYPE')
## But this still fails (argh):
## exons(txdb, columns='TXTYPE')

