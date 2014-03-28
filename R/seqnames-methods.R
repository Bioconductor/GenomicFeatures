###############################################################################
## The methods in this file are all about using the seqnames.db package in
## conjunction with Transcript.Db objects.

determineDefaultSeqnameStyle <- function(x) {
    .Deprecated("seqlevelsStyle")
    seqlevelsStyle(x)
}


## Tests:
## library(TxDb.Athaliana.BioMart.plantsmart12);
## txdb = TxDb.Athaliana.BioMart.plantsmart12; seqlevelsStyle(txdb)

## library(TxDb.Hsapiens.UCSC.hg19.knownGene);
## txdb = TxDb.Hsapiens.UCSC.hg19.knownGene; seqlevelsStyle(txdb)

