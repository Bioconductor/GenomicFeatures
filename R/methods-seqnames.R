###############################################################################
## The methods in this file are all about using the seqnames.db package in
## conjunction with Transcript.Db objects.

## Helper method to allow extraction of seqnames from TranscriptDb
setMethod("seqnames", "TranscriptDb",
          function(x){
            as.character(t(AnnotationDbi:::dbEasyQuery(
                                               AnnotationDbi:::dbConn(x),
                                               "SELECT chrom FROM chrominfo")))
          }
)


## testSeqnameStyle helper
.determineDefaultSeqnameStyle <- function(x){
  ## we will have to extract the organsim name  
  ## species <-  as.character(AnnotationDbi:::dbEasyQuery(
  ##               AnnotationDbi:::dbConn(x),
  ##               "SELECT value FROM metadata WHERE name='Genus and Species'"))
  species <- species(x)
  ## and we have to extract the seqnames (from the TxDb)
  ## seqnames <- as.character(t(AnnotationDbi:::dbEasyQuery(
  ##                            AnnotationDbi:::dbConn(x),
  ##                            "SELECT chrom FROM chrominfo")))
  
  seqnames <- seqnames(x)  ## TODO: add a helper like this.
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

## I need a method that allows me to determine what the default SeqnameStyle
## is for a TranscriptDb..
setMethod("determineDefaultSeqnameStyle", "TranscriptDb",
          function(x) .determineDefaultSeqnameStyle(x)
)


## Tests:
## library(TxDb.Athaliana.BioMart.plantsmart12); txdb = TxDb.Athaliana.BioMart.plantsmart12; determineDefaultSeqnameStyle(txdb)

## library(TxDb.Hsapiens.UCSC.hg19.knownGene); txdb = TxDb.Hsapiens.UCSC.hg19.knownGene; determineDefaultSeqnameStyle(txdb)



## seqnameStyle sets and gets the seqnameStyle string in TranscriptDb object.
## TODO: getter should be smarter.  It should check whether or not
## x$seqnameStyle is empty, and if it is, then it should set the value to the
## default using determineDefaultSeqnameStyle(txdb)
setMethod("seqnameStyle", "TranscriptDb",
          function(x) x$seqnameStyle
)

setReplaceMethod("seqnameStyle", "TranscriptDb",
    function(x, value)
    {
      species <- species(txdb)
        if (!is.null(value) && length(value==1) &&
            isSupportedSeqnamesStyle(species, value)) {
            x$seqnameStyle <- value
        }else{
          warning("That value is not a supported seqnameStyle.")
        }
        x
    }
)


## Tests:
## library(TxDb.Athaliana.BioMart.plantsmart12); txdb = TxDb.Athaliana.BioMart.plantsmart12; seqnameStyle(txdb);
## seqnameStyle(txdb) <- "UCSC" ## Should fail;
## seqnameStyle(txdb) <- "ensembl";
## seqnameStyle(txdb)

