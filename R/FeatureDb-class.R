### =========================================================================
### FeatureDb objects
### -------------------------------------------------------------------------

## This is to try and tidy up before setRefClass()
gc()

#' FeatureDb objects
#'
#' WARNING: The FeatureDb/makeFeatureDbFromUCSC/features code base is no longer
#' actively maintained and FeatureDb-related functionalities might get
#' deprecated in the near future. Please use
#' \code{\link{makeFeatureDbFromUCSC}} for a convenient way to import
#' transcript annotations from UCSC online resources into Bioconductor.
#'
#' The FeatureDb class is a generic container for storing genomic locations of
#' an arbitrary type of genomic features.
#'
#' See \code{?\link{TxDb}} for a container for storing transcript annotations.
#'
#' See \code{?\link{makeFeatureDbFromUCSC}} for a convenient way to make
#' FeatureDb objects from BioMart online resources.
#'
#'
#' @aliases FeatureDb-class class:FeatureDb FeatureDb
#' @section Methods: In the code snippets below, \code{x} is a FeatureDb
#' object.
#'
#' \describe{ \item{}{ \code{metadata(x)}: Return \code{x}'s metadata in a data
#' frame.  } }
#' @author Marc Carlson
#' @seealso \itemize{ \item The \link{TxDb} class for storing transcript
#' annotations.  \item \code{\link{makeFeatureDbFromUCSC}} for a convenient way
#' to make a FeatureDb object from UCSC online resources.  \item
#' \code{\link{saveDb}} and \code{\link{loadDb}} for saving and loading the
#' database content of a FeatureDb object.  \item \code{\link{features}} for
#' how to extract genomic features from a FeatureDb object.  }
#' @keywords methods classes
#' @examples
#'
#' fdb_file <- system.file("extdata", "FeatureDb.sqlite",
#'                         package="GenomicFeatures")
#' fdb <- loadDb(fdb_file)
#' fdb
#'
#' @exportClass FeatureDb
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

#' Extract simple features from a FeatureDb object
#'
#' WARNING: The FeatureDb/makeFeatureDbFromUCSC/features code base is no longer
#' actively maintained and FeatureDb-related functionalities might get
#' deprecated in the near future. Please use
#' \code{\link{makeFeatureDbFromUCSC}} for a convenient way to import
#' transcript annotations from UCSC online resources into Bioconductor.
#'
#' Generic function to extract genomic features from a FeatureDb object.
#'
#'
#' @aliases features features,FeatureDb-method
#' @param x A \link{FeatureDb} object.
#' @return a GRanges object
#' @author M. Carlson
#' @seealso \link{FeatureDb}
#' @examples
#'
#'   fdb <- loadDb(system.file("extdata", "FeatureDb.sqlite",
#'                                    package="GenomicFeatures"))
#'   features(fdb)
#'
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor (not exported).
###

FeatureDb <- function(conn)
{
    .FeatureDb$new(conn=conn)
}

