
#####################################################################
## Helpers to access/process the table names and columns
.getTableColMapping <- function(x){
  conn <- AnnotationDbi:::dbConn(x)
  tables <- DBI:::dbListTables(conn)
  tCols <- sapply(tables, DBI:::dbListFields, conn=conn)
  ## right up front we are getting rid of metadata and chrominfo tables...
  tCols[names(tCols)!="metadata" & names(tCols)!="chrominfo"]
}

## used to match and generate abbreviated column names
.makeColAbbreviations <- function(x){
  tCols <- .getTableColMapping(x)
  longNames <- unique(unlist(tCols,use.names=FALSE))
  abbrev <- unique(toupper((gsub("_","",unlist(tCols,use.names=FALSE)))))
  names(abbrev) <- longNames
  abbrev
}

## For when you need to get the true table names back from the abbrev's
.reverseColAbbreviations <- function(x, cnames){
  abr <- .makeColAbbreviations(x)
  names(abr)[abr %in% cnames] 
}

## used to retrieve vector of table names that go with vector of col names
.getTableNames <- function(x, cnames){
  realColNames <- .reverseColAbbreviations(x, cnames)
  tCols <- .getTableColMapping(x)
  ## Translate all names to one unique vector.
  getTabNames <- function(name, tCols){names(tCols[grep(name, tCols)])}
  tabNames <- lapply(realColNames, getTabNames, tCols)
  names(tabNames) <- realColNames
  tabNames
}

.getSimpleTableNames <- function(x, cnames){
  unique(unlist(.getTableNames(x, cnames)))
}

#####################################################################
## Helpers to generate SQL statements for select()
## g.gene_id, s.exon_rank
.makeSelectList <- function(x, cnames){
  tNames <- .getTableNames(x, cnames)
  ## For just this fun, drop all but 1st element for each tNames
  tNames <- lapply(tNames,function(x){x[1]})
  ## then continue on...  
  tabAbbrevs <- substr(unlist(tNames),1,1)
  names(tabAbbrevs) <- rep(names(tNames),elementLengths(tNames))  
  paste( paste(tabAbbrevs, ".",names(tabAbbrevs), sep=""), collapse=", ")  
}

## genes AS g, splicing AS s etc.
.makeAsList <- function(x, cnames){
  simpTNames <- .getSimpleTableNames(x, cnames)
  paste(simpTNames, "AS", substr(simpTNames,1,1), ",",collapse=" ")
}

## WHERE g._tx_id = s._tx_id etc.
.makeJoinList <- function(x, cnames){
  simpTNames <- .getSimpleTableNames(x, cnames)
  joins <- character()
  ## loop through elements of simpTNames
  for(tName in seq_len(length(simpTNames))){
    switch(EXPR = tName,
           "gene" =  joins <- c(joins, "g._tx_id = s._tx_id"),
           "transcript"  = joins <- c(joins, "t._tx_id = s._tx_id"),
           "exon"  = joins <- c(joins, "e._exon_id = s._exon_id"),
           "cds"  = joins <- c(joins, "c._cds_id = s._cds_id"))
  }
  paste(joins, collapse=" AND ")
}

.makeKeyList <- function(x, keys, keytype){
  colType <- .reverseColAbbreviations(x, keytype)
  keys <- paste(paste("'",keys,"'",sep=""),collapse=",")
  paste(colType, "IN (", keys,")")
}

.select <- function(x, keys, cols, keytype){
  ## we add TXID to cnames, which forces splicing to always be included
  ## Splicing is a almost always needed, but almost never requested.
  cnames <- c(cols, "TXID") 
  sql <- paste("SELECT",
               .makeSelectList(x, cnames),
               "FROM",
               .makeAsList(x, cnames),
               "WHERE",
               .makeJoinList(x, cnames),
               "AND",
               .makeKeyList(x, keys, keytype))
  res <- AnnotationDbi:::dbQuery(AnnotationDbi:::dbConn(x), sql)
  ## Then drop any cols that were not explicitely requested but that may have
  ## been appended to make a joind (like TXID)
  res <- res[,.reverseColAbbreviations(x,cols)]
  ## Then put the user preferred headers onto the table
  res
}


setMethod("select", "TranscriptDb",
    function(x, keys, cols, keytype) {
          .select(x, keys, cols, keytype)
        }
)




#####################################################################
## method for cols()
.cols <- function(x){
  ## cast is just to drop names here
  as.character(.makeColAbbreviations(x))  
}


setMethod("cols", "TranscriptDb",
    function(x) .cols(x)
)



#####################################################################
## method for keys()
## PHEW, There are a lot of types.  But names are not always available so...
.makeKeytypeChoiceAndGetKeys <- function(x, keytype){
    switch(EXPR = keytype,
           "GENEID" = AnnotationDbi:::dbQuery(AnnotationDbi:::dbConn(x),
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
      if (missing(keytype)) keytype <- "GENEID"
      .makeKeytypeChoiceAndGetKeys(x, keytype)
    }
)

#####################################################################
## method for keytypes()
## Seven! keytypes for TranscriptDb
setMethod("keytypes", "TranscriptDb",
    function(x) return(c("GENEID","TXID","TXNAME","EXONID","EXONNAME","CDSID","CDSNAME"))
)






##   library(TxDb.Hsapiens.UCSC.hg19.knownGene); x <- txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene; cnames = c("GENEID","TXNAME"); k = head(keys(x,keytype="GENEID")); keytype = "GENEID"; cols = c("GENEID","TXNAME")


##   keytypes(txdb)
##   head(keys(txdb, "GENEID")) 
##   head(keys(txdb, "TXID"))
##   head(keys(txdb, "EXONID"))
##   head(keys(txdb, "CDSID"))
##   head(keys(txdb, "TXNAME"))
##   head(keys(txdb, "EXONNAME"))
##   head(keys(txdb, "CDSNAME"))

##   cols(txdb)


## AnnotationDbi:::debugSQL()
##   select(txdb, head(keys(txdb, "GENEID")), cols = c("GENEID","TXNAME"), keytype="GENEID") 
