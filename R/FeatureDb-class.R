### =========================================================================
### FeatureDb objects
### -------------------------------------------------------------------------


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
    conn <- AnnotationDbi:::dbConn(x)
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Saving. 
###

setMethod("saveFeatures", "FeatureDb",
          function(x, file)
          {
            saveDb(x, file)
          }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors.  
###

## setMethod("metadata", "FeatureDb",
##     function(x) dbReadTable(conn(x), "metadata")
## )


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method for FeatureDb objects.
###

## setMethod("show", "FeatureDb",
##     function(object)
##     {
##         cat("FeatureDb object:\n")
##         metadata <- metadata(object)
##         for (i in seq_len(nrow(metadata))) {
##             cat("| ", metadata[i, "name"], ": ", metadata[i, "value"],
##                 "\n", sep="")
##         }
##     }
## )




