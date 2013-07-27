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
    seqlevels(txdb, force=TRUE) <- c(chr5 = "5")
    checkTrue(length(seqinfo(txdb))==1)
    
    ## This should work
    txdb <- restoreSeqlevels(txdb)
    seqlevels(txdb, force=TRUE) <- c(chr5 = "5", chr6="6", chr4="4")
    checkTrue(length(seqinfo(txdb))==3)
    checkIdentical(c('5','6','4'), seqlevels(txdb))
    checkTrue(seqlengths(txdb)[2] == min(seqlengths(txdb)))
    checkTrue(seqlengths(txdb)[3] == max(seqlengths(txdb)))
    
    ## And this should NOT work
    txdb <- restoreSeqlevels(txdb)
    checkException(seqlevels(txdb, force=TRUE) <- c(foo = "2"))
}


test_noChange_circ <- function(){
    txdb <- restoreSeqlevels(txdb)
    foo = seqinfo(txdb)
    foo@is_circular = rep(TRUE, 93)
    ## This should throw an exception
    checkException(seqinfo(txdb, new2old=1:93) <- foo)    
}


test_noChange_genome <- function(){
    txdb <- restoreSeqlevels(txdb)
    foo = seqinfo(txdb)
    foo@genome = rep("hg18", 93)
    ## This should throw an exception
    checkException(seqinfo(txdb, new2old=1:93) <- foo)
}


test_noChange_lengths <- function(){
    txdb <- restoreSeqlevels(txdb)
    foo = seqinfo(txdb)
    foo@seqlengths = rep(1000L, 93)
    ## This should throw an exception
    checkException(seqinfo(txdb, new2old=1:93) <- foo)
}


test_transcripts_accessors <- function(){
    txdb <- restoreSeqlevels(txdb)
    txs1 <- transcripts(txdb)
    seqlevels(txs1, force=TRUE) <- c(chr5 = "5")
    ## Then change seqlevels for txdb
    seqlevels(txdb, force=TRUE) <- c(chr5 = "5")
    txs2 <- transcripts(txdb)
    checkIdentical(txs1, txs2)

##     txdb <- restoreSeqlevels(txdb)
##     txs3 <- transcriptsBy(txdb, by="gene")
##     seqlevels(txs3, force=TRUE) <- c(chr5 = "5")
##     ## Then change seqlevels for txdb
##     seqlevels(txdb) <- c(chr5 = "5")
##     txs4 <- transcriptsBy(txdb, by="gene")
##     checkIdentical(txs3, txs4)  ## TROUBLE!!
    
}


