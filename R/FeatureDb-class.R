### =========================================================================
### FeatureDb objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A low-level accessor (not exported).
###

featuredbConn <- function(featuredb) .getConn(featuredb@envir)


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
    .valid.colnames(conn, tablename, colnames)
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
    conn <- txdbConn(x)
    c(.valid.metadata.table(conn,"FeatureDb"), 
      .valid.feature.table(conn))
}


setValidity2("FeatureDb", .valid.FeatureDb)




### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor (not exported).
###

FeatureDb <- function(conn)
{
    if (!is(conn, "SQLiteConnection"))
        stop("'conn' must be an SQLiteConnection object")
    envir <- new.env(parent=emptyenv())
    assign("conn", conn, envir=envir)
    reg.finalizer(envir, function(e) dbDisconnect(.getConn(e)))
    new("FeatureDb", envir=envir)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Saving/loading. 
###

setMethod("saveFeatures", "FeatureDb",
          function(x, file)
          {
            if (!is(x, "FeatureDb"))
              stop("'x' must be a FeatureDb object")
            if (!isSingleString(file))
              stop("'file' must be a single string")
            sqliteCopyDatabase(txdbConn(x), file)
          }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors.  
###

setMethod("metadata", "FeatureDb",
    function(x) dbReadTable(featuredbConn(x), "metadata")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method for FeatureDb objects.
###

setMethod("show", "FeatureDb",
    function(object)
    {
        cat("FeatureDb object:\n")
        metadata <- metadata(object)
        for (i in seq_len(nrow(metadata))) {
            cat("| ", metadata[i, "name"], ": ", metadata[i, "value"],
                "\n", sep="")
        }
    }
)





