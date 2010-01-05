### NOT exported (to avoid collisions with AnnotationDbi).

### Taken from AnnotationDbi (trying to avoid depending on AnnotationDbi
### for now).
.dbFileConnect <- function(dbfile)
{
    if (!file.exists(dbfile))
        stop("DB file '", dbfile, "' not found")
    #library(RSQLite)
    dbConnect(SQLite(), dbname=dbfile, cache_size=64000, synchronous=0)
}

get_dbfile <- function(libname, pkgname)
{
    filename <- paste(pkgname, ".sqlite", sep="")
    system.file("extdata", filename, package=pkgname, lib.loc=libname)
}

get_dbconn <- function(libname, pkgname)
{
    dbfile <- get_dbfile(libname, pkgname)
    .dbFileConnect(dbfile)
}

get_cached_dbfile <- function(datacache)
{
    if (!exists("dbfile", envir=datacache))
        stop("symbol \"dbfile\" not found in 'datacache'")
    get("dbfile", envir=datacache)
}

get_cached_dbconn <- function(datacache)
{
    if (!exists("dbconn", envir=datacache))
        stop("symbol \"dbconn\" not found in 'datacache'")
    get("dbconn", envir=datacache)
}

get_dbtable <- function(tablename, datacache)
{
    if (!exists(tablename, envir=datacache)) {
        dbtable <- dbReadTable(get_cached_dbconn(datacache), tablename, row.names=NULL)
        assign(tablename, dbtable, envir=datacache)
    }
    get(tablename, envir=datacache)
}

### Low-level util for setting the class of (all or some of) the columns of
### a data.frame. Could go in a low-level infrastructure package but we don't
### really have anything like that at the moment.
### Typical use:
###   x <- setDataFrameColClass(x, c(col2="integer", col5="factor"))
setDataFrameColClass <- function(x, col2class, drop.extra.cols=FALSE)
{
    if (!is.data.frame(x))
        stop("'x' must be a data.frame")
    if (!is.character(col2class) || is.null(names(col2class)))
        stop("'col2class' must be a named character vector")
    if (!all(names(col2class) %in% colnames(x)))
        stop("'col2class' has invalid names")
    if (!isTRUEorFALSE(drop.extra.cols))
        stop("'drop.extra.cols' must be TRUE or FALSE")
    y <- lapply(names(col2class),
                function(colname)
                {
                    class <- col2class[[colname]]
                    if (identical(class, "factor"))
                        as.factor(x[[colname]])
                    else
                        as(x[[colname]], class)
                }
         )
    if (drop.extra.cols) {
        names(y) <- names(col2class)
        return(as.data.frame(y, stringsAsFactors=FALSE))
    }
    x[names(col2class)] <- y
    x
}

test_GenomicFeatures <- function(dir)
{
    if (missing(dir)) {
        dir <- system.file("tests", package="GenomicFeatures")
    }
    require("RUnit", quietly=TRUE) || stop("RUnit package not found")
    suite <- defineTestSuite(name="GenomicFeatures RUnit Tests", dirs=dir,
                             testFileRegexp=".*_test\\.R$",
                             rngKind="default",
                             rngNormalKind="default")
    result <- runTestSuite(suite)
    printTextProtocol(result, showDetails=FALSE)
    if (.any_errors(result) || .any_fail(result)) {
        stop("test_GenomicFeatures FAIL")
    }
    result
}

.any_errors <- function(result) {
    any(sapply(result, function(r) r$nErr > 0))
}

.any_fail <- function(result) {
    any(sapply(result, function(r) r$nFail > 0))
}
