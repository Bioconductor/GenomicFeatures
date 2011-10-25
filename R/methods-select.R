.select <- function(x){
  sql <- paste("SELECT * FROM genes")
  AnnotationDbi:::dbQuery(AnnotationDbi:::dbConn(x), sql, 1L)
}



setMethod("select", "TranscriptDb",
    function(x, keys, cols, keytype) {
          if (missing(keytype)) keytype <- "ENTREZID"
          .select(x, keys, cols, keytype, jointype="gene_id")
        }
)


## get the list of things that can be returned.
.getTableColMapping <- function(x){
  conn <- AnnotationDbi:::dbConn(x)
  tables <- DBI:::dbListTables(conn)
  tCols <- sapply(tables, DBI:::dbListFields, conn=conn)
  tCols[names(tCols)!="metadata" & names(tCols)!="chrominfo"]
}

.cols <- function(x){
  tCols <- .getTableColMapping(x)
  unique(toupper((gsub("_","",unlist(tCols,use.names=FALSE)))))
}


setMethod("cols", "TranscriptDb",
    function(x) .cols(x)
)





## PHEW, that's a lot of types.  But names are not always available so...
## ALSO, Genes are being called a "name" to make it more intuitive for users
## who will associate "name" with "it came from someone else".  Internally we
## know that they are the only reliable name feature, but users will just
## become confused if we say gene_id when they are usually filled with names.
.makeKeytypeChoiceAndGetKeys <- function(x, keytype){
    switch(EXPR = keytype,
           "GENENAME" = AnnotationDbi:::dbQuery(AnnotationDbi:::dbConn(x),
             "SELECT DISTINCT gene_id FROM gene", 1L),
           "TXID" =  AnnotationDbi:::dbQuery(AnnotationDbi:::dbConn(x),
             "SELECT DISTINCT _tx_id FROM transcript", 1L),
           "TXNAME" =  AnnotationDbi:::dbQuery(AnnotationDbi:::dbConn(x),
             "SELECT DISTINCT tx_name FROM transcript", 1L),
           "EXONID" =  AnnotationDbi:::dbQuery(AnnotationDbi:::dbConn(x),
             "SELECT DISTINCT _exon_id FROM exon", 1L),
           "EXONNAME" =  AnnotationDbi:::dbQuery(AnnotationDbi:::dbConn(x),
             "SELECT DISTINCT exon_name FROM exon", 1L),
           "CDSID" =  AnnotationDbi:::dbQuery(AnnotationDbi:::dbConn(x),
             "SELECT DISTINCT _cds_id FROM cds", 1L),                  
           "CDSNAME" =  AnnotationDbi:::dbQuery(AnnotationDbi:::dbConn(x),
             "SELECT DISTINCT cds_name FROM cds", 1L),                  
           stop(paste(keytype, "is not a supported keytype.",
                      " Please use the keytypes",
                      "method to identify viable keytypes")))
}

## Get the list of possible keys, for a given keytype
setMethod("keys", "TranscriptDb",
    function(x, keytype){
      if (missing(keytype)) keytype <- "TXNAME"
      .makeKeytypeChoiceAndGetKeys(x, keytype)
    }
)

## Only four keytypes for TranscriptDb
setMethod("keytypes", "TranscriptDb",
    function(x) return(c("GENENAME","TXID","TXNAME","EXONID","EXONNAME","CDSID","CDSNAME"))
)






##   library(TxDb.Hsapiens.UCSC.hg19.knownGene); x <- txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

##   keytypes(txdb)
##   head(keys(txdb, "GENENAME")) 
##   head(keys(txdb, "TXID"))
##   head(keys(txdb, "EXONID"))
##   head(keys(txdb, "CDSID"))
##   head(keys(txdb, "TXNAME"))
##   head(keys(txdb, "EXONNAME"))
##   head(keys(txdb, "CDSNAME"))

##   cols(txdb)
