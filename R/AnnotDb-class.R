### =========================================================================
### AnnotDb objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A low-level accessor (not exported).
###

anndbConn <- function(anndb) .getConn(anndb@envir)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity of a AnnotDb object.
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

.valid.annot.table <- function(conn)
{
  ## Restrict column name checking to just columns that we are demanding
    colnames <- c("chrom", "strand","chromStart","chromEnd")
    msg <- .validate.colnames(conn, colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}


.valid.AnnotDb <- function(x)
{
    conn <- txdbConn(x)
    c(.valid.metadata.table(conn,"AnnotDb"), 
      .valid.annot.table(conn))
}


setValidity2("AnnotDb", .valid.AnnotDb)




### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor (not exported).
###

AnnotDb <- function(conn)
{
    if (!is(conn, "SQLiteConnection"))
        stop("'conn' must be an SQLiteConnection object")
    envir <- new.env(parent=emptyenv())
    assign("conn", conn, envir=envir)
    reg.finalizer(envir, function(e) dbDisconnect(.getConn(e)))
    new("AnnotDb", envir=envir)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Saving/loading. 
###

setMethod("saveFeatures", "AnnotDb",
          function(x, file)
          {
            if (!is(x, "AnnotDb"))
              stop("'x' must be a AnnotDb object")
            if (!isSingleString(file))
              stop("'file' must be a single string")
            sqliteCopyDatabase(txdbConn(x), file)
          }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors.  
###

setMethod("metadata", "AnnotDb",
    function(x) dbReadable(anndbConn(x), "metadata")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method for AnnotDb objects.
###

setMethod("show", "AnnotDb",
    function(object)
    {
        cat("AnnotDb object:\n")
        metadata <- metadata(object)
        for (i in seq_len(nrow(metadata))) {
            cat("| ", metadata[i, "name"], ": ", metadata[i, "value"],
                "\n", sep="")
        }
    }
)





