### makeAnnotDbFromUCSC() expects a UCSC table to have at
### least the following columns:
.UCSC_GENERICCOL2CLASS <- c(
    chrom="factor",
    strand="factor",
    chromStart="integer",
    chromEnd="integer"
)


## The following writes the contents of a generic table
.writeGenericAnnotTable <- function(conn, table, tableName, otherCols)
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
    ## TODO: switch this to the one where we use named inserts
    sqlVals <- paste("$", names(otherCols), ",", sep="")
    sqlVals[length(sqlVals)] <- sub(",","",sqlVals[length(sqlVals)])
    sql <- paste(c("INSERT INTO ",tableName,
                 " VALUES ($chrom,$strand,$chromStart,$chromEnd,",
                 sqlVals,")"), collapse="")
    dbEasyPreparedQuery(conn, sql, table)
}

## I will need a function to actually make the DB
makeAnnotDb <- function(table, tableName, otherCols)
                             ##metadata=NULL, ...)
{
    ## Create the db in a temp file.
    conn <- dbConnect(SQLite(), dbname="")
    .writeGenericAnnotTable(conn, table, tableName, otherCols)
##    .writeMetadataTable(conn, metadata)  # must come last!
    ## TODO: have to make a new object!
    AnnotDb(conn) 
}


## standard columns are chrom, chromStart, chromEnd and strand
## all others need to be specified 


makeAnnotDbFromUCSC <- function(genome="hg18",
         tablename="oreganno",
         url="http://genome.ucsc.edu/cgi-bin/",
         goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
         otherCols = c(id="character",name="character"))
{
    ## TODO: A lot of checking code below must be shared with
    ## makeTranscriptDbFromUCSC - Make sure we do all necessary checking!
    if (!isSingleString(genome))
        stop("'genome' must be a single string")
    if (!isSingleString(tablename))
        stop("'tablename' must be a single string")
    ## TODO: Modify this as we don't have our new tracks in this list (we have
    ## to make a new list)
    ## track <- supportedUCSCAnnotables()[tablename, "track"]
    track <- tablename ## temporarily just use oreganno to test
    if (is.na(track))
        stop("table \"", tablename, "\" is not supported")
    if (!isSingleString(url))
        stop("'url' must be a single string")
    if (!isSingleString(goldenPath_url))
        stop("'goldenPath_url' must be a single string")
  
    ## Create an UCSC Genome Browser session.
    ## TODO: implement supportedUCSCTables for these kinds of simple tracks
    session <- browserSession(url=url)
    genome(session) <- genome
    track_tables <- tableNames(ucscTableQuery(session, track=track))
    if (!(tablename %in% track_tables))
        stop("GenomicFeatures internal error: ", tablename, " table doesn't ",
             "exist or is not associated with ", track, " track. ",
             "Thank you for reporting this to the GenomicFeatures maintainer ",
             "or to the Bioconductor mailing list, and sorry for the ",
             "inconvenience.")

    ## Download the transcript table.
    message("Download the ", tablename, " table ... ", appendLF=FALSE)
    query <- ucscTableQuery(session, track, table=tablename,
                            names=NULL)
    ucsc_table <- getTable(query)
    
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
    
    message("Make the AnnoDb object ... ", appendLF=FALSE)
    makeAnnotDb(table=ucsc_table, tableName=tablename,
                ## metadata=metadata,
                otherCols)


}



## Test Code:
## library(GenomicFeatures)
## foo = makeAnnotDbFromUCSC()
## saveFeatures(foo, "annoDb.sqlite")


