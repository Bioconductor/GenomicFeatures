## Very simple extractor to just return what is in our Feature.Db objects

## business end of things

.extractDataCols <- function(conn, tableName){
  SQL <- paste("SELECT * FROM ",tableName,";",sep="")
  ans <- dbEasyQuery(conn, SQL)
}


.extractFeaturesAsGRanges <- function(db){
  ## 1st figure out what table is not the metadata table.
  conn <- AnnotationDbi:::dbConn(db) ## featuredbConn(db)
  tableNames <- dbListTables(conn)
  tableName <- tableNames[!tableNames %in% "metadata"]

  ## Then learn what the columns are in that table and assign to otherCols
  colNames <- dbListFields(conn, tableName)
  colNames <- colNames[!colNames %in%
                      c("chrom","strand","chromStart","chromEnd")]

  ## Extract the data
  df <- .extractDataCols(conn, tableName)
  
  ## Make & return the Object
  ranges <- IRanges(start= unlist(df["chromStart"]),
                    end = unlist(df["chromEnd"]))
  ans <- GRanges(seqnames = unlist(df["chrom"]),
                 ranges = ranges,
                 strand = unlist(df["strand"]))
  otherVals <- as(df[colNames], "DataFrame")
  values(ans) <- otherVals
  metadata(ans)[[1]] <- DataFrame(metadata(db))
  ans
}

## method:
setMethod("features", "FeatureDb",
          function(x){.extractFeaturesAsGRanges(x)}
)

## test code:
## library(GenomicFeatures)
## fdb <- loadDb("FeatureDb.sqlite")
## features(fdb)
