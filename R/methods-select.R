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
## Helpers for join-type selection

## this just takes the 1 letter abrevs and makes them into a sorted string
## that can be used as a key below
.encodeSortedTableKey <- function(sTNames){
  prefSort <- c("g","t","s","e","c")  
  res <- sTNames[match(prefSort, sTNames)]
  paste(res[!is.na(res)], collapse="")
}
.makeTableKey <- function(x,cnames){
  sTNames <- substr(.getSimpleTableNames(x, cnames),1,1)
  .encodeSortedTableKey(sTNames)    
}

## for unlikely table combos 
.missingTableInterpolator <- function(tName){
  tName <- switch(EXPR = tName,
                  "se" = "tse",
                  "sc" = "tsc",
                  "te" = "tse",
                  "tc" = "tsc",
                  "ge" = "gtse",
                  "gc" = "gtsc",
                  "gs" = "gts",
                  "sec" = "tsec",
                  "gsec" = "gtsec",
                  "gtce" = "gtsec",
                  "gte" = "gtse",
                  tName)
  tName
}

## real joins for likely combos
.tableJoinSelector <- function(tName){
  ## if its not one of these, then it needs to become one
  tName <- .missingTableInterpolator(tName)
  gt <- paste("SELECT * FROM transcript LEFT JOIN gene ",
                 "ON (transcript._tx_id = gene._tx_id) ")
  gts <- paste("SELECT * FROM transcript LEFT JOIN gene ",
                 "ON (transcript._tx_id = gene._tx_id) INNER JOIN splicing ",
                 "ON (transcript._tx_id = splicing._tx_id) ")
  gtse <- paste("SELECT * FROM transcript LEFT JOIN gene ",
                 "ON (transcript._tx_id = gene._tx_id) INNER JOIN splicing ",
                 "ON (transcript._tx_id = splicing._tx_id) ",
                 "INNER JOIN exon ON (splicing._exon_id = exon._exon_id) ")
  gtsc <- paste("SELECT * FROM transcript LEFT JOIN gene ",
                 "ON (transcript._tx_id = gene._tx_id) INNER JOIN splicing ",
                 "ON (transcript._tx_id = splicing._tx_id) ",
                 "LEFT JOIN cds ON (splicing._cds_id = cds._cds_id) ")
  gtsec <- paste("SELECT * FROM transcript LEFT JOIN gene ",
                 "ON (transcript._tx_id = gene._tx_id) INNER JOIN splicing ",
                 "ON (transcript._tx_id = splicing._tx_id) ",
                 "INNER JOIN exon ON (splicing._exon_id = exon._exon_id) ",
                 "LEFT JOIN cds ON (splicing._cds_id = cds._cds_id) ")
  ts <- paste("SELECT * FROM transcript INNER JOIN splicing ",
              "ON (transcript._tx_id = splicing._tx_id) ")
  tse <- paste("SELECT * FROM transcript INNER JOIN splicing ",
               "ON (transcript._tx_id = splicing._tx_id) ",
               "INNER JOIN exon ON (splicing._exon_id = exon._exon_id) ")
  tsc <- paste("SELECT * FROM transcript INNER JOIN splicing ",
               "ON (transcript._tx_id = splicing._tx_id) ",
               "LEFT JOIN cds ON (splicing._cds_id = cds._cds_id) ")
  tsec <- paste("SELECT * FROM transcript INNER JOIN splicing ",
                "ON (transcript._tx_id = splicing._tx_id) ",
                "INNER JOIN exon ON (splicing._exon_id = exon._exon_id) ",
                "LEFT JOIN cds ON (splicing._cds_id = cds._cds_id) ")
  
  sql <- switch(EXPR = tName,
                "g" = "SELECT * FROM gene as g",
                "t" = "SELECT * FROM transcript as t",
                "s" = "SELECT * FROM splicing as s",
                "e" = "SELECT * FROM exon as e",
                "c" = "SELECT * FROM cds as c",
                "gt" = gt,
                "gts" = gts,
                "gtse" = gtse,
                "gtsc" = gtsc,
                "gtsec" = gtsec,
                "tse" = tse,
                "tsc" = tsc,
                "tsec" = tsec,
                stop(paste("No query for this combination of tables.",
                           "Please add",tName,"to the interpolator")))
  sql
}


#####################################################################
## Helpers to generate SQL statements for select()

## g.gene_id, s.exon_rank
## For some cols, they will occur in more than one table within the join,
## in that case, we just grab the 1st one.
.makeSelectList <- function(x, cnames, abbrev=TRUE){
    tNames <- .getTableNames(x, cnames)
    ## Here is where we only grab the 1st one...
    tNames <- lapply(tNames,function(x){x[1]})
    ## then continue on...  
    tabAbbrevs <- substr(unlist(tNames),1,1)
    names(tabAbbrevs) <- rep(names(tNames),elementLengths(tNames))    
  if(abbrev==TRUE){
    paste( paste(tabAbbrevs, ".",names(tabAbbrevs), sep=""), collapse=", ")
  }else{
    paste(names(tabAbbrevs), collapse=", ")
  }
}

## genes AS g, splicing AS s etc.
.makeAsList <- function(x, cnames){
  simpTNames <- .getSimpleTableNames(x, cnames)
  paste(simpTNames, "AS", substr(simpTNames,1,1), collapse=", ")
}

## WHERE g._tx_id = s._tx_id etc.
.makeJoinSQL <- function(x, cnames){
  tKey <- .makeTableKey(x,cnames)
  .tableJoinSelector(tKey)  
}


.makeKeyList <- function(x, keys, keytype, abbrev=TRUE){
  #colType <- .reverseColAbbreviations(x, keytype)
  colType <- .makeSelectList(x, keytype, abbrev)
  keys <- paste(paste("'",keys,"'",sep=""),collapse=",")
  paste(colType, "IN (", keys,")")
}

.select <- function(x, keys, cols, keytype){
  ## 1st we check the keytype to see if it is valid:
  if(is.na(keys(x, keytype)[1]) & length(keys(x, keytype))==1){ 
    stop(paste("There do not appear to be any keys",
               "for the keytype you have specified."))
    }
  ## we used to add TXID to cnames, which forces splicing to always be included
  ## Splicing is a almost always needed, but almost never requested.
  ## cnames <- unique(c(cols, "TXID", keytype))
  ## 
  cnames <- unique(c(cols, keytype))
  tKey <- .makeTableKey(x,cnames)
#message(paste("keytype generated:",tKey))
  sql <- paste("SELECT DISTINCT",
               .makeSelectList(x, cnames, abbrev=FALSE),
               "FROM (",
               .makeJoinSQL(x, cnames),
               ") WHERE",
               .makeKeyList(x, keys, keytype, abbrev=FALSE))
  res <- AnnotationDbi:::dbQuery(AnnotationDbi:::dbConn(x), sql)
  ## Then sort rows and cols and drop the filtered rows etc. using .resort
  ## from AnnoationDbi
  joinType <- .reverseColAbbreviations(x, keytype)
  if(dim(res)[1]>0){
    res <- AnnotationDbi:::.resort(res, keys, joinType)
  }
  ## Then drop any cols that were not explicitely requested but that may have
  ## been appended to make a joind (like TXID)
  res <- res[,.reverseColAbbreviations(x,cols),drop=FALSE]
##TODO: implement the above (equiv.) for the AnnotationDbi .select()
  
  ## Then I need to filter out rows of NAs
  res <- res[!apply(is.na(res),1,all),,drop=FALSE]
  ## always reset rownames after removing rows
  rownames(res) <- NULL
  
  ## Then put the user preferred headers onto the table
  fcNames <- .makeColAbbreviations(x)
  colnames(res) <- fcNames[match(colnames(res), names(fcNames))]
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
  res <- .makeColAbbreviations(x)
  names(res) <- NULL
  res
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






##   library(TxDb.Hsapiens.UCSC.hg19.knownGene); x <- txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene;

##   cols = c("GENEID"); keys = head(keys(x, "GENEID")); foo = select(x, keys, cols = cols, keytype="GENEID");head(foo)

##   cols = c("GENEID","TXID"); keys = head(keys(x, "GENEID")); foo = select(x, keys, cols = cols, keytype="GENEID");head(foo)

##   cols = c("GENEID","TXID", "EXONRANK"); keys = head(keys(x, "GENEID")); foo = select(x, keys, cols = cols, keytype="GENEID");head(foo)



##   cols = c("GENEID","TXID", "EXONRANK","CDSID"); keys = head(keys(x, "GENEID")); foo = select(x, keys, cols = cols, keytype="GENEID");head(foo)


##   cols = c("GENEID","TXID", "EXONRANK", "EXONID"); keys = head(keys(x, "GENEID")); foo = select(x, keys, cols = cols, keytype="GENEID");head(foo)


## SUPER SLOW:
##   cols = c("GENEID","TXID", "EXONRANK", "EXONID", "CDSID"); keys = head(keys(x, "GENEID")); foo = select(x, keys, cols = cols, keytype="GENEID");head(foo)

## TODO: WHY is THIS a gtsec? Where is gene coming from?  (I fear it is from
## TXID...)
##   cols = c("TXID", "EXONRANK", "EXONID", "CDSID"); keys = head(keys(x, "TXID")); foo = select(x, keys, cols = cols, keytype="TXID");head(foo)




##   cols = c("EXONRANK", "EXONID", "CDSID"); keys = head(keys(x, "EXONID")); foo = select(x, keys, cols = cols, keytype="EXONID");head(foo)












##  cnames = c("GENEID","TXNAME"); k = head(keys(x,keytype="GENEID")); keytype = "GENEID"; cols = c("GENEID","TXNAME")


##   keytypes(txdb)
##   head(keys(txdb, "GENEID")) 
##   head(keys(txdb, "TXID"))
##   head(keys(txdb, "EXONID"))
##   head(keys(txdb, "CDSID"))
##   head(keys(txdb, "TXNAME"))
##   head(keys(txdb, "EXONNAME"))
##   head(keys(txdb, "CDSNAME"))

##   cols(txdb)


##   AnnotationDbi:::debugSQL()  
##   k = c("foo",head(keys(txdb, "GENEID"))); foo = select(txdb, k, cols = c("GENEID","TXNAME"), keytype="GENEID"); head(foo)


## more testing:
##   cols = cols(x); k = head(keys(x,keytype="GENEID")); foo = select(txdb, k, cols = cols, keytype="GENEID");head(foo)





## A SQL bug seems to happen if I ask for a bit of everything (every table)
## and provide a key that is not an entrez gene ID

## cols = c("GENEID","TXNAME", "CDSID", "EXONSTART"); k = head(keys(txdb, "TXNAME")); foo = select(txdb, k, cols = cols, keytype="TXNAME"); head(foo)

## what is odd also is that the SQL (in a SQL session) does not seem to
## complete...


## It's not just exonstart stuff because this works:
##  cols = c("GENEID", "EXONSTART","TXNAME"); k = head(keys(txdb, "TXNAME")); foo = select(txdb, k, cols = cols, keytype="TXNAME"); head(foo)


## While this fails
##  cols = c("GENEID", "CDSID","TXNAME"); k = head(keys(txdb, "TXNAME")); foo = select(txdb, k, cols = cols, keytype="TXNAME"); head(foo)




















###############################################
############ Cleared up
###############################################
## Also this gives me a strange answer (one that is not really filtered quite right) - There may be a problem with AnnotationDbi:::.resort() - OR maybe I can just put a "DISTINCT" into the query? - and maybe I should really do both???

## cols = c("GENEID","TXNAME", "TXID"); k = head(keys(txdb, "TXID"));foo = select(txdb, k, cols = cols, keytype="TXID"); head(foo)

## Issue (used to) exist with above where the txid is NOT actually being
## filtered
## right.  The answer that comes back looks like:
##   TXID GENEID     TXNAME
## 1    1   <NA>       <NA>
## 2    2   <NA>       <NA>
## 3    3   <NA>       <NA>
## 4    4 653635 uc009vis.2
## 5    4 653635 uc009vis.2
## 6    4 653635 uc009vis.2
## ## but should look like:
##   TXID GENEID     TXNAME
## 1    1   <NA>       <NA>
## 2    2   <NA>       <NA>
## 3    3   <NA>       <NA>
## 4    4 653635 uc009vis.2
## 5    5 653635 uc009vjc.1
## 6    6 653635 uc009vjd.2


## DONE: add checks to pre-warn when a keytype has no values for it's key
## cols = c("GENEID","TXNAME", "TXID", "CDSNAME"); k = head(keys(txdb, "CDSNAME"));foo = select(txdb, k, cols = cols, keytype="CDSNAME"); head(foo)




## The following breaks because I end up with a jointype of _cds_id, but with
## NO COLLUMNS that match to it..  So my SQL code has failed to retrieve
## everything I needed...

## IOW, the res has content from the DB but not the right amount?
## cols = c("GENEID","TXNAME", "TXID"); k = head(keys(txdb, "CDSID"));foo = select(txdb, k, cols = cols, keytype="CDSID"); head(foo)


## Another example of this:
## The following wasn't working because the key column must always be part of
## the query. (fixed)

##  debug(AnnotationDbi:::.resort); AnnotationDbi:::debugSQL();
##  cols = c("GENEID","TXNAME", "TXID", "CDSNAME"); k = head(keys(txdb, "CDSID"));foo = select(txdb, k, cols = cols, keytype="CDSID"); head(foo)

## TODO: Drop rows that contain only NAs. - This is a good cleanup generalyl, and also if you drop the column that you sorted on then you may have added a bunch of NA rows for where there was no data, and now those are just silly. - DONE






#### So what is the pattern?
#### prefix (add later)
## SELECT * FROM

#### OPTIONAL: swap in as clause code for "gene AS g" etc.
#### OPTIONAL EFFICIENCY: replace * from core clauses with string from .makeSelectList() (but after filtering out unwanted stuff of course) so that we only grab stuff that was asked for from the beginning.
#### OPTIONAL EFFICIENCY???: Can I make it smarter so that it starts with the table that has the selected keys and then builds out from there?  No that can't work, because of the outer join requirements?? - definitely want to do this one "later" since even if it can work it will be difficult to implement gracefully.

#### create 8 core-combos
#### and select appropriate one based on who is in .getSimpleTableNames()
## (SELECT * FROM gene AS g LEFT JOIN transcript AS t ON
## g._tx_id=t._tx_id
## UNION
## SELECT * FROM transcript AS t LEFT JOIN gene AS g ON
## g._tx_id=t._tx_id )

#### suffix (add later, once we have our cores)
## WHERE gene_id IN ('10772','22947')
## limit 10;


