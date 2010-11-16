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

## helper function to correct for UCSC data having off by one start info.
.adjustchromStarts <- function(table){
  chromStart <- as.integer(table[["chromStart"]])
  table <- table[,!(colnames(table) %in% "chromStart")]
  chromStart <- chromStart + 1L
  return(cbind(table,chromStart))
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
.writeGenericFeatureTable <- function(conn, table, tableName, columns)
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
    sql2 <- paste("  ", names(columns), " TEXT,\n")
    sql <- c(sql1, sql2,")")
    ## remove final comma
    sql[length(sql)-1] <- sub(",","",sql[length(sql)-1])
    dbEasyQuery(conn, paste(sql, collapse=""))
    ## Fill the  table.
    sqlVals <- paste("$", names(columns), ",", sep="")
    sqlVals[length(sqlVals)] <- sub(",","",sqlVals[length(sqlVals)])
    sql <- paste(c("INSERT INTO ",tableName,
                 " VALUES ($chrom,$strand,$chromStart,$chromEnd,",
                 sqlVals,")"), collapse="")
    dbEasyPreparedQuery(conn, sql, table)
}




## TODO: expose these and document them.  Do I want the arg name to be "track"
## or tablename for UCSCFeatureDbTrackSchemas() ???
supportedUCSCFeatureDbTracks <- function(genome="hg19")
{
  ##TODO: fill out black list of tracks that we cannot support here.
  ##unsupported <- c("ruler")
  unsupported <- c("fakeTrackName")
  session <- browserSession()
  genome(session) <- genome
  query <- ucscTableQuery(session)
  trackNames(query)[!(trackNames(query) %in% unsupported)]
}


UCSCFeatureDbTrackSchemas <- function(genome="hg19", track="stsMap")
{
  session <- browserSession()
  genome(session) <- genome  
  query <- ucscTableQuery(session, track)  
  res <- ucscSchema(query)
  ## now for the tricky part: converting from MYSQL to R...  There is no good
  ## way to extract the "R" type information from the data.frame since it
  ## appears that they are all treated as "character" information.  And there
  ## is no good way to convert from MySQL to R since that is handled elsewhere
  ## at the C-level and baked into the code that extracts the data.  So here
  ## we will just assume (as michael does for his $example slot) that he has
  ## done something reasonable in making all these things to be "character"
  ## (and realistically, this is probably fine for what we are doing here)
  ## So for now, just fake up the fake track info. from rtracklayer...
  types <- unlist(lapply(res@listData$example, class))
  names <- res@listData$field  
  result <- types
  names(result) <- names
  result
}




## I will need a function to actually make the DB
makeFeatureDb <- function(table, tableName, columns, metadata=NULL, ...)
{
    ## Create the db in a temp file.
    conn <- dbConnect(SQLite(), dbname="")
    .writeGenericFeatureTable(conn, table, tableName, columns)
    .writeMetadataFeatureTable(conn, metadata, tableName)  # must come last!
    FeatureDb(conn) 
}



## standard columns are chrom, chromStart, chromEnd and strand
## all others need to be specified 
makeFeatureDbFromUCSC <- function(genome="hg18",
         tablename="oreganno",
         columns = UCSCFeatureDbTrackSchemas(genome, tablename),
         url="http://genome.ucsc.edu/cgi-bin/",
         goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
{
    if (!isSingleString(genome))
        stop("'genome' must be a single string")
    if (!isSingleString(tablename))
        stop("'tablename' must be a single string")
    
    ## Extract the tablename provided from the list of viable tracks.   
    trks <- supportedUCSCFeatureDbTracks(genome)
    track <- trks[trks %in% tablename]      
    if (length(track)==0)
        stop("table \"", tablename, "\" is not supported")

    ## Check the column names
    if(length(names(columns)) != length(unique(names(columns))))
      stop("The default field names are not unique for this table.")
    ## Once we know the columns names are unique, we remove the default ones.
    columns <- columns[!(names(columns) %in% names(.UCSC_GENERICCOL2CLASS))] 

    ## Check other arguments
    if (!isSingleString(url))
        stop("'url' must be a single string")
    if (!isSingleString(goldenPath_url))
        stop("'goldenPath_url' must be a single string")

    ## Create a UCSC Genome Browser session.
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
    ucsc_table <- .adjustchromStarts(ucsc_table)
    
    ## check that we have at least the 5 columns of data
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
    .UCSC_GENERICCOL2CLASS = c(.UCSC_GENERICCOL2CLASS, columns)
    ucsc_table <- setDataFrameColClass(ucsc_table ,.UCSC_GENERICCOL2CLASS,
                                     drop.extra.cols=TRUE)    
    
    ## Compile some of the metadata
    metadata <- .prepareUCSCFeatureMetadata(genome, tablename)
    
    message("Make the AnnoDb object ... ", appendLF=FALSE)
    makeFeatureDb(table=ucsc_table, tableName=tablename,
                metadata=metadata,
                columns)
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
##                  columns=c(name="character",gieStain="character"))




## The following is a great example
## library(GenomicFeatures)
## foo = makeFeatureDbFromUCSC("hg18","cytoBand")


## TODO: see if I can figure out why some of these tracks don't work...

## like "vegaGeneComposite" and "encodeYaleAffyRNATars"

