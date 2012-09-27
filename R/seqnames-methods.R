###############################################################################
## The methods in this file are all about using the seqnames.db package in
## conjunction with Transcript.Db objects.

## testSeqnameStyle helper
.determineDefaultSeqnameStyle <- function(x){
  ## we will have to extract the organsim name  
  species <- species(x)
  ## and we have to extract the seqnames (from the TxDb)
  seqnames <- seqlevels(x)
  ## Then for all the possible styles, pick one that works and return it...
  styles <- supportedSeqnameStyles()[[sub(" ","_",species)]]
  styleIdx <- testSeqnames(styles=styles,seqnames=seqnames,species=species)
  styles <- styles[styleIdx]
  if(length(styles) > 0){
    return(styles)[1] ## just return the 1st "true" one
    ## If there were multiple matches the distinction would not be relevant
  }else{
    warning("No styles in the seqnames.db database match what is in this TranscriptDb database.  It might be worthwhile to add a the new style to that database.")
    return(NA)
  }
}

setGeneric("determineDefaultSeqnameStyle",
           function(x) standardGeneric("determineDefaultSeqnameStyle"))

setMethod("determineDefaultSeqnameStyle", "TranscriptDb",
          function(x) .determineDefaultSeqnameStyle(x)
)

## Tests:
## library(TxDb.Athaliana.BioMart.plantsmart12); txdb = TxDb.Athaliana.BioMart.plantsmart12; determineDefaultSeqnameStyle(txdb)

## library(TxDb.Hsapiens.UCSC.hg19.knownGene); txdb = TxDb.Hsapiens.UCSC.hg19.knownGene; determineDefaultSeqnameStyle(txdb)

.setseqnameStyle <- function(x, value){
  species <- species(txdb)
  if (!is.null(value) && length(value==1) &&
      isSupportedSeqnamesStyle(style=value, species=species)) {
    x$seqnameStyle <- value
  }else{
    warning("That value is not a supported seqnameStyle.")
  }
  x
}

#setReplaceMethod("seqnameStyle", "TranscriptDb",
#          function(x,value) .setseqnameStyle(x,value)
#)

## Tests:
## library(TxDb.Athaliana.BioMart.plantsmart12); txdb = TxDb.Athaliana.BioMart.plantsmart12; seqnameStyle(txdb);
## seqnameStyle(txdb) <- "UCSC" ## Should fail;
## seqnameStyle(txdb) <- "NCBI";
## seqnameStyle(txdb)

