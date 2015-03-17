## 1st we need to write tests for all the helper functions then we need more
## tests for select (generally)
## Why test the lower level helpers?  Because that way I will get a failure
## point right at the location where the trouble occurs (high resolution for
## trouble detection)._
require("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require("RUnit")
  
test_getTableColMapping <- function(){
  res <- GenomicFeatures:::.getTableColMapping(txdb)
  exp <- list(cds=c("_cds_id","cds_name","cds_chrom","cds_strand","cds_start",
                "cds_end"),
              exon=c("_exon_id","exon_name","exon_chrom","exon_strand",
                "exon_start","exon_end"),
              gene=c("gene_id","_tx_id"),
              splicing=c("_tx_id","exon_rank","_exon_id","_cds_id"),
              transcript=c("_tx_id","tx_name","tx_type","tx_chrom","tx_strand",
                "tx_start","tx_end"))
  checkIdentical(res, exp)
}

test_makeColAbbreviations <- function(){
  res <- GenomicFeatures:::.makeColAbbreviations(txdb)
  checkTrue(res[["_cds_id"]]=="CDSID")
  res2 <- GenomicFeatures:::.getTableColMapping(txdb)
  checkTrue(length(res)==21, length(unique(unlist(res2))))
}

test_reverseColAbbreviations <- function(){
  cnames <- GenomicFeatures:::.makeColAbbreviations(txdb)
  res <- GenomicFeatures:::.reverseColAbbreviations(txdb, cnames)
  checkTrue(names(cnames)[[1]]==res[[1]])
  checkTrue(length(res) == length(cnames))
}

test_getTableNames <- function(){
  cnames <- GenomicFeatures:::.makeColAbbreviations(txdb)
  res <- GenomicFeatures:::.getTableNames(txdb, cnames)
  ## Let's check the ones that fail more easily
  checkTrue(length(res[["_tx_id"]])==3)
  checkTrue(length(res[["_exon_id"]])==2)
  checkTrue(length(res[["_cds_id"]])==2)
}

test_getSimpleTableNames <- function(){
  cnames <- GenomicFeatures:::.makeColAbbreviations(txdb)
  res <- GenomicFeatures:::.getSimpleTableNames(txdb, cnames)
  exp <- c("cds","splicing","exon","gene","transcript")
  checkIdentical(res, exp)
}

test_encodeSortedTableKey <- function(){
  sTNames <- c("s", "e", "t")
  res <- GenomicFeatures:::.encodeSortedTableKey(sTNames)
  checkIdentical(res,"tse")
}

test_makeTableKey <- function(){
  cnames <- GenomicFeatures:::.makeColAbbreviations(txdb)
  res <- GenomicFeatures:::.makeTableKey(txdb, cnames)
  checkIdentical(res,"gtsec")
}

test_missingTableInterpolator <- function(){
  tName <- "x"
  res <- GenomicFeatures:::.missingTableInterpolator(tName)
  checkIdentical(res,"x")
  tName <- "se"
  res <- GenomicFeatures:::.missingTableInterpolator(tName)
  checkIdentical(res,"tse")
}

test_tableJoinSelector <- function(){
  tName <- "t"
  res <- GenomicFeatures:::.tableJoinSelector(tName)
  checkIdentical(res,"transcript")
  tName <- "gt"
  res <- GenomicFeatures:::.tableJoinSelector(tName)
  checkIdentical(res, "(SELECT * FROM transcript LEFT JOIN gene  ON (transcript._tx_id = gene._tx_id) )")
  tName <- "FOOBAZZLE"
  checkException(GenomicFeatures:::.tableJoinSelector(tName))
}

test_makeSelectList <- function(){
  cnames <- GenomicFeatures:::.makeColAbbreviations(txdb)[c("_cds_id","_tx_id")]
  res <- GenomicFeatures:::.makeSelectList(txdb, cnames)
  exp <- "c._cds_id, g._tx_id" ## 2nd one will be a "g."
  checkIdentical(res, exp)
}

test_makeAsList <- function(){
  cnames <- GenomicFeatures:::.makeColAbbreviations(txdb)[c("_cds_id","_tx_id")]
  res <- GenomicFeatures:::.makeAsList(txdb, cnames)
  exp <- "cds AS c, splicing AS s, gene AS g, transcript AS t"
  checkIdentical(res, exp)
}

test_makeJoinSQL <- function(){
  cnames <- GenomicFeatures:::.makeColAbbreviations(txdb)[c("_cds_id","_tx_id")]
  res <- GenomicFeatures:::.makeJoinSQL(txdb, cnames)
  exp <- "(SELECT * FROM transcript LEFT JOIN gene  ON (transcript._tx_id = gene._tx_id) INNER JOIN splicing  ON (transcript._tx_id = splicing._tx_id)  LEFT JOIN cds ON (splicing._cds_id = cds._cds_id) )"
  checkIdentical(res, exp)
}

test_makeKeyList <- function(){
  ks <- 1:6
  kt <- "TXID"
  res <- GenomicFeatures:::.makeKeyList(txdb, keys=ks, keytype=kt)
  exp <- "g._tx_id IN ( '1','2','3','4','5','6' )"
  checkIdentical(res, exp)
}


test_keys <- function(){
  checkException(keys(txdb, keytype="CDSCHROM"))
}

test_keys_advancedArgs <- function(){
    k1 <- keys(txdb, keytype="TXNAME")
    checkTrue("uc001aaa.3" %in% k1)
    
    k2 <- keys(txdb, keytype="TXNAME", pattern=".2$")
    checkTrue("uc001aaq.2" %in% k2)
    checkTrue(!("uc001aaa.3" %in% k2))
    checkTrue(length(k2) < length(k1))

    l1 <- length(keys(txdb, keytype="TXID", column="GENEID"))
    l2 <- length(keys(txdb, keytype="TXID"))
    checkTrue(l1 < l2)
    
    k3 <- head(keys(txdb, keytype="GENEID", pattern=".2$",
                    column="TXNAME", fuzzy=TRUE))
    res <- suppressWarnings( select(txdb, k3, columns=c("GENEID","TXNAME"),
                                   keytype="GENEID"))
    checkTrue(any(grepl(".2$",res$TXNAME)))
}



test_select <- function(){
  keys = head(keys(txdb, "GENEID"))
  cols = c("GENEID")
  res <- select(txdb, keys, cols, keytype="GENEID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("GENEID"), colnames(res))
  checkTrue(length(res$GENEID)==length(keys))
  checkIdentical(res$GENEID, keys)

  keys = head(keys(txdb, "TXID"))
  cols = c("TXID")
  res <- select(txdb, keys, cols, keytype="TXID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("TXID"), colnames(res))
  checkTrue(length(res$TXID)==length(keys))
  checkIdentical(res$TXID, keys)
 
  keys = head(keys(txdb, "GENEID"))
  cols = c("GENEID","TXID")
  res <- select(txdb, keys, cols, keytype="GENEID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("GENEID","TXID"), colnames(res))
  checkTrue(length(unique(res$GENEID))==length(keys))

  keys = head(keys(txdb, "GENEID"))
  cols = c("GENEID","TXID", "EXONRANK")
  res <- select(txdb, keys, cols, keytype="GENEID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("GENEID","TXID","EXONRANK"), colnames(res))
  checkTrue(length(unique(res$GENEID))==length(keys))

  keys = head(keys(txdb, "GENEID"))
  cols = c("GENEID","TXID", "EXONRANK","CDSID")
  res <- select(txdb, keys, cols, keytype="GENEID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("GENEID","CDSID","TXID","EXONRANK"), colnames(res))
  checkTrue(length(unique(res$GENEID))==length(keys))

  ## It's really cosmetic but: should the order of the final data.frame match
  ## the order of the cols?
  ## I think so, except that we may add a col for keys (even if not requested)
  ## if added, such a col should be in front.
  
  keys = head(keys(txdb, "GENEID"))
  cols = c("GENEID","TXID", "EXONRANK", "EXONID")
  res <- select(txdb, keys, cols, keytype="GENEID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("GENEID","EXONID","TXID","EXONRANK"), colnames(res))
  checkTrue(length(unique(res$GENEID))==length(keys))

  keys = head(keys(txdb, "GENEID"))
  cols = c("GENEID","TXID", "EXONRANK", "EXONID", "CDSID")
  res <- select(txdb, keys, cols, keytype="GENEID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("GENEID","CDSID","EXONID","TXID","EXONRANK"), colnames(res))
  checkTrue(length(unique(res$GENEID))==length(keys))
  
  keys = head(keys(txdb, "TXID"))
  cols = c("TXID", "EXONRANK", "EXONID", "CDSID")
  res <- select(txdb, keys, cols, keytype="TXID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("TXID","CDSID","EXONID","EXONRANK"), colnames(res))
  checkTrue(length(unique(res$TXID))==length(keys))

  keys = head(keys(txdb, "EXONID"))
  cols = c("EXONRANK", "EXONID", "CDSID")
  res <- select(txdb, keys, cols, keytype="EXONID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("EXONID","CDSID","EXONRANK"), colnames(res))
  checkTrue(length(unique(res$EXONID))==length(keys))
  
  keys = head(keys(txdb, "TXNAME"))
  cols = c("GENEID","TXNAME", "CDSID", "EXONSTART")
  res <- select(txdb, keys, cols, keytype="TXNAME")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("TXNAME","CDSID","EXONSTART","GENEID"), colnames(res))  
  checkTrue(length(unique(res$TXNAME))==length(keys))

  
  keys = head(keys(txdb, "TXNAME"))
  cols = c("GENEID", "EXONSTART","TXNAME")
  res <- select(txdb, keys, cols, keytype="TXNAME")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("TXNAME","EXONSTART","GENEID"), colnames(res))
  checkTrue(length(unique(res$TXNAME))==length(keys))
    
  
  keys = head(keys(txdb, "TXNAME"))
  cols = c("GENEID", "CDSID","TXNAME")
  res <- select(txdb, keys, cols, keytype="TXNAME")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("TXNAME","CDSID","GENEID"), colnames(res))
  checkTrue(length(unique(res$TXNAME))==length(keys))
    
  
  keys = head(keys(txdb, "TXID"))
  cols = c("GENEID","TXNAME", "TXID")
  res <- select(txdb, keys, cols, keytype="TXID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("TXID","GENEID","TXNAME"), colnames(res))
  checkTrue(length(unique(res$TXID))==length(keys))
  ## For this particular case, we want to make sure that the TXNAMES are not
  ## being copied (there should be one unique one for each ID in this range)
  checkTrue(length(unique(res$TXNAME)) == length(res$TXNAME))
  
  keys = head(keys(txdb, "CDSNAME"))
  cols = c("GENEID","TXNAME", "TXID", "CDSNAME")
  checkException(select(txdb, keys, cols, keytype="CDSNAME"), silent=TRUE)
  
  keys = head(keys(txdb, "CDSID"))
  cols = c("GENEID","TXNAME", "TXID", "CDSNAME")
  res <- select(txdb, keys, cols, keytype="CDSID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols)+1) ## this is one where we ADD an extra!
  checkIdentical(c("CDSID","CDSNAME","GENEID","TXID","TXNAME"), colnames(res))
  checkTrue(length(unique(res$CDSID))==length(keys))

  
  ## stress test (this used to take way too long)
  keys = head(keys(txdb, "GENEID"))
  cols = c("GENEID","CDSSTART")
  res <- select(txdb, keys, cols, keytype="GENEID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("GENEID","CDSSTART"), colnames(res))
    
}

test_select_isActiveSeq <- function(){
  
  ## set isActiveSeq to only watch chr1
 txdb <- restoreSeqlevels(txdb)  ## This is to reset things (safety measure)
 isActiveSeq(txdb)[seqlevels(txdb)] <- FALSE
 isActiveSeq(txdb) <- c("chr1"=TRUE)  
  
  ## then use select
  keys <- head(keys(txdb, "GENEID"))
  cols <- c("GENEID","CDSSTART", "CDSCHROM")
  res <- select(txdb, keys, columns = cols, keytype="GENEID")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("GENEID","CDSCHROM","CDSSTART"), colnames(res))
  uniqChrs <- unique(res$CDSCHROM)[!is.na(unique(res$CDSCHROM))]
  checkIdentical(c("chr1"),uniqChrs)

  ## keys must contain keys that match to more than one thing
  keys <- c(head(keys(txdb,keytype="TXNAME")),
            tail(keys(txdb,keytype="TXNAME")))
  cols <- c("TXNAME","TXCHROM","TXSTRAND")
  res <- select(txdb, keys, columns = cols, keytype="TXNAME")
  checkTrue(dim(res)[1]>0)
  checkTrue(dim(res)[2]==length(cols))
  checkIdentical(c("TXNAME","TXCHROM","TXSTRAND"), colnames(res))
  uniqChrs <- unique(res$TXCHROM)[!is.na(unique(res$TXCHROM))]
  checkIdentical(c("chr1"),uniqChrs)  
}



