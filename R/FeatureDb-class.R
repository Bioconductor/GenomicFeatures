### =========================================================================
### FeatureDb objects
### -------------------------------------------------------------------------


## This is to try and tidy up before setRefClass()
gc()

.FeatureDb <-
    setRefClass("FeatureDb", contains="AnnotationDb")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A low-level accessor (not exported).
###

## featuredbConn <- function(featuredb) featuredb$conn


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity of a FeatureDb object.
###

.validate.colnames <- function(conn, colnames)
{
    ## even though we don't know the name of the table, take advantage of the
    ## fact that there are only two tables and one of them is always called
    ## "metadata"
    tablenames <- dbListTables(conn)
    tablename <- tablenames[!tablenames %in% "metadata"]
    AnnotationDbi:::.valid.colnames(conn, tablename, colnames)
}

.valid.feature.table <- function(conn)
{
  ## Restrict column name checking to just columns that we are demanding
    colnames <- c("chrom", "strand","chromStart","chromEnd")
    msg <- .validate.colnames(conn, colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}


.valid.FeatureDb <- function(x)
{
    conn <- dbconn(x)
    c(AnnotationDbi:::.valid.metadata.table(conn, "Db type",
                                            "FeatureDb"),
      .valid.feature.table(conn))
}


setValidity2("FeatureDb", .valid.FeatureDb)




### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor (not exported).
###

FeatureDb <- function(conn)
{
    .FeatureDb$new(conn=conn)
}

