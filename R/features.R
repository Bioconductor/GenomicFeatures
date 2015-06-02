## Very simple extractor to just return what is in our Feature.Db objects

## business end of things

.extractDataCols <- function(conn, tableName){
  SQL <- paste0("SELECT * FROM ", tableName, ";")
  dbEasyQuery(conn, SQL)
}


.extractFeaturesAsGRanges <- function(db)
{
    ## 1st figure out what table is not the metadata table.
    conn <- dbconn(db) ## featuredbconn(db)
    tableNames <- dbListTables(conn)
    tableName <- tableNames[!tableNames %in% "metadata"]

    ## Then learn what the columns are in that table and assign to otherCols
    colNames <- dbListFields(conn, tableName)
    reserved <- c("name", "chrom", "strand", "chromStart", "chromEnd")
    colNames <- colNames[!colNames %in% reserved]

    ## Extract the data
    df <- .extractDataCols(conn, tableName)
    
    ## Make & return the Object
    md <- metadata(db)
    genome <- md[md$name == "Genome", 'value']
    if (is.null(genome))
        genome <- NA_character_
    ans <-
        GRanges(seqnames = df$chrom,
                ranges = IRanges(df$chromStart, df$chromEnd, names=df$name),
                strand = sub('\\.', '*', df$strand),
                df[colNames],
                seqinfo = Seqinfo(unique(df$chrom), genome=genome))
    metadata(ans)[[1]] <- DataFrame(metadata(db))
    ans
}

setGeneric("features", signature="x",
           function(x) standardGeneric("features"))

setMethod("features", "FeatureDb",
          function(x) .extractFeaturesAsGRanges(x))

## test code:
## library(GenomicFeatures)
## fdb <- loadDb("FeatureDb.sqlite")
## features(fdb)
