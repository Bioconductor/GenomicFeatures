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
  chromStart <- chromStart + 1L
  table$chromStart <- chromStart
  table
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
    organism <- lookup_organism_by_UCSC_genome(genome)
    if (is.na(taxonomyId)) {
        taxonomyId <- GenomeInfoDb:::lookup_tax_id_by_organism(organism)
    } else {
        GenomeInfoDb:::check_tax_id(taxonomyId)
    }
    
    metadata <- data.frame(
        name=c("Data source", "Genome", "UCSC Table", "Organism",
               "Taxonomy ID", "Resource URL"),
        value=c("UCSC", genome, tablename, organism,
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
.writeGenericFeatureTable <- function(conn, data, tableName, columns)
{
    data <- unique(data)
    ## for now just drop lines that don't have values for chromStart
    ## we may need to be more stringent
    data <- S4Vectors:::extract_data_frame_rows(data, !is.na(data$chromStart))
    ## Create the table.
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
    dbExecute(conn, paste(sql, collapse=""))
    ## Fill the table.
    insert_data_into_table(conn, tableName, data)
}


## helper function to ID tables that rtracklayer won't process.
checkTable <- function(query){
  query@table %in% rtracklayer:::tableNames(query)
}

## helper to check a track
isGoodTrack <- function(track, session){
  query <- ucscTableQuery(session)
  tracks <- trackNames(query)
  track %in% names(tracks)
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
  tables <- ucscTables(genome, track)
  if(length(tables)){
    tables
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
  query <- ucscTableQuery(session, table=tablename)
  res <- ucscSchema(query)
  ## now for the tricky part: converting from MYSQL to R...  There is no good
  ## way to extract the "R" type information from the data.frame since it
  ## appears that they are all treated as "character" information.
  sqlTypes <- res$SQL.type
  Rtypes <- map_SQLtypes_to_Rtypes(sqlTypes)
  names <- res$field
  names(Rtypes) <- names
  Rtypes
}

## Convert SQL types to R types by creating a
## dummy SQL table and reading it's type information
map_SQLtypes_to_Rtypes <- function(SQLtypes)
{
    stopifnot(is.character(SQLtypes))
    SQLtypes <- gsub(" *unsigned", "", SQLtypes)
    col_defs <- sprintf("col%d %s", seq_along(SQLtypes), SQLtypes)
    sql <- sprintf("CREATE TABLE dummy (%s)", paste0(col_defs, collapse=", "))
    conn <- dbConnect(SQLite())
    on.exit(dbDisconnect(conn))
    dbExecute(conn, sql)
    dummy <- dbReadTable(conn, "dummy")
    col_types <- vapply(dummy, function(col) class(col)[[1L]], character(1))
    # edge case
    col_types[col_types == "blob"] <- "character"
    setNames(col_types, SQLtypes)
}

## I will need a function to actually make the DB
makeFeatureDb <- function(data, tableName, columns, metadata=NULL, ...)
{
    ## Create the db in a temp file.
    conn <- dbConnect(SQLite(), dbname="")
    .writeGenericFeatureTable(conn, data, tableName, columns)
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
         goldenPath.url=getOption("UCSC.goldenPath.url"),
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
    if (!isSingleString(goldenPath.url))
        stop("'goldenPath.url' must be a single string")

    ## Create a UCSC Genome Browser session.
    session <- browserSession(url=url)
    genome(session) <- genome
    track_tables <- ucscTables(genome, track)
    if (!(tablename %in% track_tables))
        stop("GenomicFeatures internal error: ", tablename, " table doesn't ",
             "exist or is not associated with ", track, " track. ",
             "Thank you for reporting this to the GenomicFeatures maintainer ",
             "or to the Bioconductor mailing list, and sorry for the ",
             "inconvenience.")
    
    ## Download the data table.
    message("Download the ", tablename, " table ... ", appendLF=FALSE)
    query <- ucscTableQuery(session, table=tablename)
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
    ## ensure that the table columns conform to expectations
    ucsc_table <- ucsc_table[,names(.UCSC_GENERICCOL2CLASS)]
    
    ## Compile some of the metadata
    metadata <- .prepareUCSCFeatureMetadata(genome, tablename, taxonomyId)
    
    message("Make the AnnoDb object ... ")
    makeFeatureDb(data=ucsc_table, tableName=tablename,
                  metadata=metadata,
                  columns)
}


