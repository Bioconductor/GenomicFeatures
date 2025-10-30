### =========================================================================
### Extractors for features in other databases
### -------------------------------------------------------------------------
###
### This is for extractors that do NOT point to the TxDb proper.
### Such extractors can point to other databases OR they can
### point to other FeatureDbs within the same package.

.syncSeqlevel <- function(txdb, ans){
  isActSeq <- .isActiveSeq(txdb)
  n2oNames <- levels(seqnames(ans))
  n2o <- match(seqnames(seqinfo(txdb)), n2oNames)
  seqinfo(ans, new2old=n2o) <- seqinfo(txdb)
  seqlevels(ans, pruning.mode="coarse") <- names(isActSeq)[isActSeq]
  ans
}

## main function
.tRNAs <- function(x) {
  S4Vectors:::load_package_gracefully("txdbmaker", "when calling ",
                                      "tRNAs() on a TxDb object")
  fdbpkg <- "FDb.UCSC.tRNAs"
  fdbenv <- loadNamespace(fdbpkg)
  ## get the current package name
  pkgName <- txdbmaker::makePackageName(x)
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

