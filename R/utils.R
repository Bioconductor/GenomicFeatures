### =========================================================================
### Miscellaneous low-level utils
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DB related.
###
### NOT exported (to avoid collisions with AnnotationDbi).
###

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Data frame related.
###

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

### Acts like an SQL *inner* join.
### 'x' must be a data frame. 'name2val' must be a named atomic vector or
### a named factor. 'join_colname' must be the name of the col in 'x' whose
### values are matched against 'names(name2val)'. 'vals_colname' must be the
### name of the col that will be populated with the appropriate 'name2val'
### vals and bound to 'x'.
### Note that this acts like an SQL *inner* join, not a *left* join, i.e.
### rows in 'x' that can't be mapped to a value in 'name2val' are dropped.
joinDataFrameWithName2Val <- function(x, join_colname, name2val, vals_colname)
{
    if (!is.data.frame(x))
        stop("'x' must be a data.frame")
    if (!isSingleString(join_colname) || is.null(x[[join_colname]]))
        stop("'join_colname' must be a valid colname for 'x'")
    if (!is.vector(name2val) && !is.factor(name2val))
        stop("'name2val' must be a vector (or factor)")
    if (!is.atomic(name2val) || is.null(names(name2val)))
        stop("'name2val' must be atomic and have names")
    if (!isSingleString(vals_colname) || !is.null(x[[vals_colname]]))
        stop("invalid 'vals_colname'")
    tmp <- split(as.vector(name2val), names(name2val))
    ## as.vector() is required below just because 'x[[join_colname]]' could
    ## be a factor (subsetting by a factor is equivalent to subsetting by
    ## an integer vector but this is not what we want here).
    tmp <- tmp[as.vector(x[[join_colname]])]
    x <- x[rep.int(seq_len(nrow(x)), elementLengths(tmp)), ]
    row.names(x) <- NULL
    if (nrow(x) == 0L)
        vals <- name2val[FALSE]
    else if (is.factor(name2val))
        vals <- factor(unname(unlist(tmp)), levels=levels(name2val))
    else
        vals <- unname(unlist(tmp))
    x[[vals_colname]] <- vals
    x
}

### Returns the vector of ids such that 'unique(x)[ids, ]' is identical
### to 'x' (in the same way that 'levels(f)[f]' is identical to
### 'as.vector(f)' when 'f' is a character factor).
### This unambiguously defines 'ids'. In particular, it's not Locale
### specific, despite the fact that the current implementation uses a
### sorting approach.
makeIdsForUniqueDataFrameRows <- function(x)
{
    if (!is.data.frame(x))
        stop("'x' must be a data.frame")
    x_order <- do.call(order, x)
    x_dups <- duplicated(x)
    ## First we make "provisory" ids. Those ids *are* Locale specific.
    prov_ids <- integer(nrow(x))
    prov_ids[x_order] <- cumsum(!x_dups[x_order])
    ## Convert the "provisory" ids into the final ids. The final ids are
    ## *not* Locale specific anymore.
    as.integer(factor(prov_ids, levels=unique(prov_ids)))
}

