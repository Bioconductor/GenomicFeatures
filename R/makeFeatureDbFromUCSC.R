### makeFeatureDbFromUCSC() expects a UCSC table to have at
### least the following columns:
.UCSC_GENERICCOL2CLASS <- c(
    chrom="factor",
    strand="factor",
    chromStart="integer",
    chromEnd="integer"
)

## helper function to add missing strand information to table (if missing)
.addMissingStrandCols <- function(table){
  if(!"strand" %in% colnames(table)){
    strand <- rep("*", dim(table)[1])
    return(cbind(table,strand))
  }else{
    return(table)
  }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### prepare and write out the DB table contents 
###

.prepareUCSCFeatureMetadata <- function(genome, tablename)
{
    message("Prepare the 'metadata' data frame ... ",
            appendLF=FALSE)
    metadata <- data.frame(
        name=c("Data source", "Genome", "UCSC Table"),
        value=c("UCSC", genome, tablename)
    )
    message("OK")
    metadata
}

### Some of the metadata is added later.
.writeMetadataFeatureTable <- function(conn, metadata, tableName)
{
    data_nrow <- dbEasyQuery(conn, paste("SELECT COUNT(*) FROM ",tableName,
                                         collapse=""))[[1L]]    
    thispkg_version <- installed.packages()['GenomicFeatures', 'Version']
    rsqlite_version <- installed.packages()['RSQLite', 'Version']
    mat <- matrix(c(
        DB_TYPE_NAME, "FeatureDb",
        "data_nrow", data_nrow,
        "Db created by", "GenomicFeatures package from Bioconductor",
        "Creation time", svn.time(),
        "GenomicFeatures version at creation time", thispkg_version,
        "RSQLite version at creation time", rsqlite_version,
        "DBSCHEMAVERSION", DB_SCHEMA_VERSION),
        ncol=2, byrow=TRUE
    )
    colnames(mat) <- c("name", "value")
    metadata <- rbind(data.frame(name=mat[ , "name"], value=mat[ , "value"],
                                 stringsAsFactors=FALSE),
                      metadata)
    dbWriteTable(conn, "metadata", metadata, row.names=FALSE)
}

## The following writes the data contents of our generic table
.writeGenericFeatureTable <- function(conn, table, tableName, otherCols)
{
    table <- unique(table)
    ## for now just drop lines that don't have values for chromStart
    ## we may need to be more stringent
    table <- table[!is.na(table$chromStart), ]
    ## Create the 'tableName' table.
    sql1 <- c("CREATE TABLE ",tableName," (\n",
              "  chrom TEXT NOT NULL,\n",
              "  strand TEXT NOT NULL,\n",
              "  chromStart INTEGER NOT NULL,\n",
              "  chromEnd INTEGER NOT NULL,\n") ## always a comma = never done
    ## Add remaining rows (b/c there will ALWAYS be at least one "other" field)
    sql2 <- paste("  ", names(otherCols), " TEXT,\n")
    sql <- c(sql1, sql2,")")
    ## remove final comma
    sql[length(sql)-1] <- sub(",","",sql[length(sql)-1])
    dbEasyQuery(conn, paste(sql, collapse=""))
    ## Fill the  table.
    sqlVals <- paste("$", names(otherCols), ",", sep="")
    sqlVals[length(sqlVals)] <- sub(",","",sqlVals[length(sqlVals)])
    sql <- paste(c("INSERT INTO ",tableName,
                 " VALUES ($chrom,$strand,$chromStart,$chromEnd,",
                 sqlVals,")"), collapse="")
    dbEasyPreparedQuery(conn, sql, table)
}

## TODO: explore using rtracklayer to discover these? 
## TODO: add code to put in * when there is no strand information.
.SUPPORTED_FEATUREDB_UCSC_TABLES <- c(
  ## tablename (unique key)    track             subtrack
  ## "",                 "",     NA,
  ## "",                 "",     NA,
  "stsMap",                 "STS Markers",     NA,
  "cytoBand",                 "Chromosome Band",     NA,
  "oreganno",                 "ORegAnno",     NA
)

supportedUCSCFeatureDbTracks <- function()
{
    mat <- matrix(.SUPPORTED_FEATUREDB_UCSC_TABLES, ncol=3, byrow=TRUE)
    colnames(mat) <- c("tablename", "track", "subtrack")
    data.frame(track=mat[ , "track"], subtrack=mat[ , "subtrack"],
               row.names=mat[ , "tablename"],
               stringsAsFactors=FALSE)
}



## I will need a function to actually make the DB
makeFeatureDb <- function(table, tableName, otherCols, metadata=NULL, ...)
{
    ## Create the db in a temp file.
    conn <- dbConnect(SQLite(), dbname="")
    .writeGenericFeatureTable(conn, table, tableName, otherCols)
    .writeMetadataFeatureTable(conn, metadata, tableName)  # must come last!
    FeatureDb(conn) 
}


## standard columns are chrom, chromStart, chromEnd and strand
## all others need to be specified 

makeFeatureDbFromUCSC <- function(genome="hg18",
         tablename="oreganno",
         url="http://genome.ucsc.edu/cgi-bin/",
         goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
         otherCols = c(id="character",name="character"))
{
    if (!isSingleString(genome))
        stop("'genome' must be a single string")
    if (!isSingleString(tablename))
        stop("'tablename' must be a single string")
    ## TODO: implement supportedUCSCTables for these kinds of simple tracks
    ## Modify this as we don't have our new tracks in this list (we have
    ## to make a new list)
    
    track <- supportedUCSCFeatureDbTracks()[tablename, "track"]    
    ## track <- tablename ## temporarily just use oreganno to test
    
    if (is.na(track))
        stop("table \"", tablename, "\" is not supported")
    if (!isSingleString(url))
        stop("'url' must be a single string")
    if (!isSingleString(goldenPath_url))
        stop("'goldenPath_url' must be a single string")
  
    ## Create an UCSC Genome Browser session.
    session <- browserSession(url=url)
    genome(session) <- genome
    track_tables <- tableNames(ucscTableQuery(session, track=track))
    if (!(tablename %in% track_tables))
        stop("GenomicFeatures internal error: ", tablename, " table doesn't ",
             "exist or is not associated with ", track, " track. ",
             "Thank you for reporting this to the GenomicFeatures maintainer ",
             "or to the Bioconductor mailing list, and sorry for the ",
             "inconvenience.")

    ## Download the data table.
    message("Download the ", tablename, " table ... ", appendLF=FALSE)
    query <- ucscTableQuery(session, track, table=tablename,
                            names=NULL)
    ucsc_table <- getTable(query)

    ## check that we have strand info, and if not, add some in
    ucsc_table <- .addMissingStrandCols(ucsc_table)
    
    ## check that we have at least the 5 columns
    if (ncol(ucsc_table) < length(.UCSC_GENERICCOL2CLASS)+1)
        stop("GenomicFeatures internal error: ", tablename, " table doesn't ",
             "exist, was corrupted during download, or doesn't contain ",
             "sufficient information. ",
             "Thank you for reporting this to the GenomicFeatures maintainer ",
             "or to the Bioconductor mailing list, and sorry for the ",
             "inconvenience.")
    message("OK")
    
    ## then make our table, but remember, we have to add new columns to our
    ## base table type 1st:
    .UCSC_GENERICCOL2CLASS = c(.UCSC_GENERICCOL2CLASS,otherCols)
    ucsc_table <- setDataFrameColClass(ucsc_table ,.UCSC_GENERICCOL2CLASS,
                                     drop.extra.cols=TRUE)    

    ## Compile some of the metadata
    metadata <- .prepareUCSCFeatureMetadata(genome, tablename)
    
    message("Make the AnnoDb object ... ", appendLF=FALSE)
    makeFeatureDb(table=ucsc_table, tableName=tablename,
                metadata=metadata,
                otherCols)


}



## Test Code:
## library(GenomicFeatures)
## foo = makeFeatureDbFromUCSC()
## saveFeatures(foo, "FeatureDb.sqlite")

## library(GenomicFeatures)
## foo = loadFeatures("FeatureDb.sqlite")


## Code to just discover the tracks (and columns in ea.)



## library(GenomicFeatures)
## foo = makeFeatureDbFromUCSC("hg18","cytoBand",
##                  otherCols=c(name="character",gieStain="character"))
