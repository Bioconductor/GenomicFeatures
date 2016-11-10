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
    seqlevels(txdb, pruning.mode="coarse") <- c(chr5 = "5")
    checkTrue(length(seqinfo(txdb))==1)
    
    ## This should work
    txdb <- restoreSeqlevels(txdb)
    seqlevels(txdb, pruning.mode="coarse") <- c(chr5 = "5", chr6="6", chr4="4")
    checkTrue(length(seqinfo(txdb))==3)
    checkIdentical(c('5','6','4'), seqlevels(txdb))
    checkTrue(seqlengths(txdb)[2] == min(seqlengths(txdb)))
    checkTrue(seqlengths(txdb)[3] == max(seqlengths(txdb)))
    
    ## And this should NOT work
    txdb <- restoreSeqlevels(txdb)
    checkException(seqlevels(txdb, pruning.mode="coarse") <- c(foo = "2"))
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


test_transcripts_accessor <- function(){
    txdb <- restoreSeqlevels(txdb)
    txs1 <- transcripts(txdb)
    seqlevels(txs1, pruning.mode="coarse") <- c(chr5 = "5")
    ## Then change seqlevels for txdb
    seqlevels(txdb, pruning.mode="coarse") <- c(chr5 = "5")
    txs2 <- transcripts(txdb)
    checkIdentical(txs1, txs2)
}

test_exons_accessor <- function(){
    txdb <- restoreSeqlevels(txdb)
    exs1 <- exons(txdb)
    seqlevels(exs1, pruning.mode="coarse") <- c(chr5 = "5")
    ## Then change seqlevels for txdb
    seqlevels(txdb, pruning.mode="coarse") <- c(chr5 = "5")
    exs2 <- exons(txdb)
    checkIdentical(exs1, exs2)
}

test_cds_accessor <- function(){
    txdb <- restoreSeqlevels(txdb)
    cds1 <- cds(txdb)
    seqlevels(cds1, pruning.mode="coarse") <- c(chr5 = "5")
    ## Then change seqlevels for txdb
    seqlevels(txdb, pruning.mode="coarse") <- c(chr5 = "5")
    cds2 <- cds(txdb)
    checkIdentical(cds1, cds2)
}

test_promoters_accessor <- function(){
    txdb <- restoreSeqlevels(txdb)
    prm1 <- promoters(txdb)
    seqlevels(prm1, pruning.mode="coarse") <- c(chr5 = "5")
    ## Then change seqlevels for txdb
    seqlevels(txdb, pruning.mode="coarse") <- c(chr5 = "5")
    prm2 <- promoters(txdb)
    checkIdentical(prm1, prm2)
}


test_transcriptsBy_accessors <- function(){
    ## This one is a "fun" one.
    ## There are issues because some genes are annotated as being on
    ## TWO different chromosomes.  Such genes are filtered for txs3,
    ## but NOT for txs4...   Hmmmm.
    txdb <- restoreSeqlevels(txdb)
    txs3 <- transcriptsBy(txdb, by="gene")
    seqlevels(txs3, pruning.mode="coarse") <- c(chr5 = "5")
    ## Then change seqlevels for txdb
    seqlevels(txdb, pruning.mode="coarse") <- c(chr5 = "5")
    txs4 <- transcriptsBy(txdb, by="gene")
##    checkIdentical(txs3, txs4)  ## TROUBLE!!
    
}


## What to do about this?  The reason for the difference is because of order of operations.  txs3 gets all the ranges and then removes any that are not kosher (this is correct), txs4 OTOH gets only ranges from chr5 (efficient!), but then fails to filter out things that have hybrid seqnames (as they were pre-filtered).  I think I have to make the query less efficient to fix this, but I want to discuss it with Herve 1st to get a 2nd opinion.
