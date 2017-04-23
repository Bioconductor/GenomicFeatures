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

## helper function to re-assign column names as required
## if the value has been set and there is precisely ONE col that matches the 
## new val, then go ahead and rename it to the chrom, chromStart or chromEnd.
.renamColsHelper <- function(table, col, newString){
          tmpColNames <- colnames(table)
          tmpColNames[colnames(table)==col] <- newString
          colnames(table) <- tmpColNames  
          table
}

.checkAndRenamCols <- function(table, chromCol, chromStartCol, chromEndCol){
    if(!is.null(chromCol) && table(colnames(table)==chromCol)[["TRUE"]]==1) {
          table <-.renamColsHelper(table, chromCol, "chrom")
    }
    if(!is.null(chromStartCol) 
        && table(colnames(table)==chromStartCol)[["TRUE"]]==1){
          table <- .renamColsHelper(table, chromStartCol, "chromStart")
    }
    if(!is.null(chromEndCol) 
        && table(colnames(table)==chromEndCol)[["TRUE"]]==1){
          table <-.renamColsHelper(table, chromEndCol, "chromEnd")
    }
  table
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### prepare and write out the DB table contents 
###

.prepareUCSCFeatureMetadata <- function(genome, tablename, taxonomyId)
{
    message("Prepare the 'metadata' data frame ... ",
            appendLF=FALSE)
    if(is.na(taxonomyId)){
        taxonomyId <- GenomeInfoDb:::.taxonomyId(UCSCGenomeToOrganism(genome))
    }else{
        GenomeInfoDb:::.checkForAValidTaxonomyId(taxonomyId)
    }
    
    metadata <- data.frame(
        name=c("Data source", "Genome", "UCSC Table", "Organism",
               "Taxonomy ID", "Resource URL"),
        value=c("UCSC", genome, tablename, UCSCGenomeToOrganism(genome),
                taxonomyId, "http://genome.ucsc.edu/")
    )
    message("OK")
    metadata
}

### Some of the metadata is added later.
.writeMetadataFeatureTable <- function(conn, metadata, tableName) {
    data_nrow <- dbEasyQuery(conn, paste("SELECT COUNT(*) FROM ",tableName,
                                         collapse=""))[[1L]]
    thispkg_version <- packageDescription("GenomicFeatures")$Version
    rsqlite_version <- packageDescription("RSQLite")$Version
    mat <- matrix(c(
        DB_TYPE_NAME, "FeatureDb",
        "Supporting package", "GenomicFeatures",
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
    table <- S4Vectors:::extract_data_frame_rows(table,
                                                 !is.na(table$chromStart))
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
    sqlVals <- paste0("$", names(columns), ",")
    sqlVals[length(sqlVals)] <- sub(",","",sqlVals[length(sqlVals)])
    sql <- paste(c("INSERT INTO ",tableName,
                 " VALUES ($chrom,$strand,$chromStart,$chromEnd,",
                 sqlVals,")"), collapse="")
    dbEasyPreparedQuery(conn, sql, table)
}


## helper function to ID tables that rtracklayer won't process.
checkTable <- function(query){
  "primaryTable" %in% rtracklayer:::ucscTableOutputs(query)
}

## helper to check a track
isGoodTrack <- function(track, session){
  query <- ucscTableQuery(session, track=track)
  checkTable(query)
}

## helper to detect and generate a list of "legitimate" tracks
makeWhiteList <- function(session, trx){
  sapply(trx, isGoodTrack, session)
}

## Discovery for supported Tracks
supportedUCSCFeatureDbTracks <- function(genome)
{
  session <- browserSession()
  genome(session) <- genome
  query <- ucscTableQuery(session)
  trx <- trackNames(ucscTableQuery(session))
  supported <- makeWhiteList(session, trx)
  trx[supported]
}

## Discover table names available in Tracks
supportedUCSCFeatureDbTables <- function(genome, track)
{
  session <- browserSession()
  genome(session) <- genome
  query <- ucscTableQuery(session, track=track)
  if(checkTable(query)){
    tableNames(query)	
  }else{
    stop("The track provided does not contain tables that are available in tabular form.")
  }
}

## Discover the schema information (field names and potentially someday the
## type information) for a track and table combo.
UCSCFeatureDbTableSchema <- function(genome,
                                     track,
                                     tablename)
{
  session <- browserSession()
  genome(session) <- genome
  ## Check that the track is available for this genome   
  if(!isGoodTrack(track, session))
    stop("track \"", track, "\" is not supported")
  ## Check that the tablename is available for this genome
  tbls <- supportedUCSCFeatureDbTables(genome, track)
  tbl <- tbls[tbls %in% tablename]
  if (length(tbl)==0)
    stop("table \"", tablename, "\" is not supported")
  
  ## then make a query
  query <- ucscTableQuery(session,
                          track=track,
                          table=tablename)  
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
  ## maybe someday i will be able to get more complete information from
  ## rtracklayer, but for now, character is ok
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
makeFeatureDbFromUCSC <- function(genome,
         track,
         tablename,
         columns = UCSCFeatureDbTableSchema(genome, track, tablename),
         url="http://genome.ucsc.edu/cgi-bin/",
         goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
         chromCol=NULL,
         chromStartCol=NULL,
         chromEndCol=NULL,
         taxonomyId=NA)
{
    if (!isSingleString(genome))
        stop("'genome' must be a single string")
    if (!isSingleString(track))
        stop("'track' must be a single string")
    if (!isSingleString(tablename))
        stop("'tablename' must be a single string")
    
    ## Check the column names
    if(length(names(columns)) != length(unique(names(columns))))
      stop("The default field names are not unique for this table.")
    ## Once we know the columns names are unique, we remove the default ones.
    columns <- columns[!(names(columns) %in% names(.UCSC_GENERICCOL2CLASS))] 
    ## also have to remove any columns that are to be re-assigned!
    if(!is.null(chromCol) || !is.null(chromStartCol) || !is.null(chromEndCol)){
      ## if I concatenate to a vector, the NULL values will be MIA = perfect
      altCols <- c(chromCol,chromStartCol,chromEndCol) 
      columns <- columns[!(names(columns) %in% altCols )]
    }

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
    query <- ucscTableQuery(session, track, table=tablename)
    ucsc_table <- getTable(query)
    
    ## check that we have strand info, and if not, add some in
    ucsc_table <- .addMissingStrandCols(ucsc_table)
    ## TODO: do any required substitutions (required column renames)
    ucsc_table <- .checkAndRenamCols(ucsc_table,
                                     chromCol,
                                     chromStartCol,
                                     chromEndCol)
    
    
    ## check that we have at least the 5 columns of data
    if(ncol(ucsc_table) < length(.UCSC_GENERICCOL2CLASS)+1)
        stop("GenomicFeatures internal error: ", tablename, " table doesn't ",
             "exist, was corrupted during download, or doesn't contain ",
             "sufficient information. ",
             "Thank you for reporting this to the GenomicFeatures maintainer ",
             "or to the Bioconductor mailing list, and sorry for the ",
             "inconvenience.")
    message("OK")
    ## also check that our data CONTAINS the column names we need it to


    message("Checking that required Columns are present ... ")
    if(length(intersect(colnames(ucsc_table),names(.UCSC_GENERICCOL2CLASS)))<4) 
        ## That means that some required cols are missing!
        stop("GenomicFeatures internal error: ", tablename, " table doesn't ",
             "contain a 'chrom', 'chromStart', or 'chromEnd' column and no ",
             "reasonable substitute has been designated via the 'chromCol'",
             "'chromStartCol' or 'chromEndCol' arguments.  ",
             "If this is not possible, please report that fact to the ",
             "GenomicFeatures maintainer or to the Bioconductor mailing list. ",
             " Thank you.")
    message("OK")


    ## Then adjust the chromStarts column (corrects for unusual UCSC counting)
    ucsc_table <- .adjustchromStarts(ucsc_table)


    ## then make our table, but remember, we have to add new columns to our
    ## base table type 1st:
    .UCSC_GENERICCOL2CLASS = c(.UCSC_GENERICCOL2CLASS, columns)
    ucsc_table <- setDataFrameColClass(ucsc_table ,.UCSC_GENERICCOL2CLASS,
                                     drop.extra.cols=TRUE)    
    
    ## Compile some of the metadata
    metadata <- .prepareUCSCFeatureMetadata(genome, tablename, taxonomyId)
    
    message("Make the AnnoDb object ... ")
    makeFeatureDb(table=ucsc_table, tableName=tablename,
                metadata=metadata,
                columns)
}


