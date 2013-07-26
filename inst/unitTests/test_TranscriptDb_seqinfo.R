library(TxDb.Hsapiens.UCSC.hg19.knownGene);
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene

test_rename_seqlevels <- function(){
    seqlevels(txdb) <- as.character(1:length(seqlevels(txdb)))
    checkIdentical(as.character(1:length(seqlevels(txdb))),
                   seqlevels(txdb))
}

test_restrict_seqlevels <- function(){
    ## This should work
    txdb <- restoreSeqlevels(txdb)    
    seqlevels(txdb) <- c(chr5 = "5")
    checkTrue(length(seqinfo(txdb))==1)
    
    ## This should work
    txdb <- restoreSeqlevels(txdb)
    seqlevels(txdb) <- c(chr5 = "5", chr6="6", chr4="4")
    checkTrue(length(seqinfo(txdb))==3)
    checkIdentical(c('5','6','4'), seqlevels(txdb))
    checkTrue(seqlengths(txdb)[2] == min(seqlengths(txdb)))
    checkTrue(seqlengths(txdb)[3] == max(seqlengths(txdb)))
    
    ## And this should NOT work
    txdb <- restoreSeqlevels(txdb)
    checkException(seqlevels(txdb) <- c(foo = "2"))
}


## test_noChange_circ <- function(){
##     txdb <- restoreSeqlevels(txdb)
##     foo = seqinfo(txdb)
##     foo@is_circular = rep(TRUE, 93)
##     ## This should throw an exception
##     checkException(seqinfo(txdb, new2old=1:93) <- foo)    
## }


## test_noChange_circ <- function(){
##     txdb <- restoreSeqlevels(txdb)
##     foo = seqinfo(txdb)
##     foo@genome = rep("hg18", 93)
##     ## This should throw an exception
##     checkException(seqinfo(txdb, new2old=1:93) <- foo)
## }

## test_noChange_circ <- function(){
##     txdb <- restoreSeqlevels(txdb)
##     foo = seqinfo(txdb)
##     foo@seqlengths = rep(1000L, 93)
##     ## This should throw an exception
##     checkException(seqinfo(txdb, new2old=1:93) <- foo)
## }
