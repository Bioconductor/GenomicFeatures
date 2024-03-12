### =========================================================================
### Miscellaneous low-level utils
### -------------------------------------------------------------------------
###
### Unless stated otherwise, nothing in this file is exported.
###


### TODO: Move this to S4Vectors (or BiocBaseUtils).
load_package_gracefully <- function(package, ...)
{
    if (!requireNamespace(package, quietly=TRUE))
        stop("Could not load package ", package, ". Is it installed?\n\n  ",
             wmsg("Note that ", ..., " requires the ", package, " package. ",
                  "Please install it with:"),
             "\n\n    BiocManager::install(\"", package, "\")")
}

call_fun_in_txdbmaker <- function(fun, ...)
{
    load_package_gracefully("txdbmaker", "Starting with BioC 3.19, ",
                            "calling ", fun, "()")
    msg <- c(fun, "() has moved to the txdbmaker package. Please ",
             "call txdbmaker::", fun, "() to get rid of this warning.")
    warning(wmsg(msg))
    FUN <- base::get(fun, envir=asNamespace("txdbmaker"), inherits=FALSE)
    do.call(FUN, list(...))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DB related.
###
### Most of this stuff was copy/pasted from AnnotationDbi (trying to avoid
### depending on AnnotationDbi for now).
###

### Environment for storing run-time objects
RTobjs <- new.env(hash=TRUE, parent=emptyenv())

assign("debugSQL", FALSE, envir=RTobjs)

debugSQL <- function()
{
    debugSQL <- !get("debugSQL", envir=RTobjs)
    assign("debugSQL", debugSQL, envir=RTobjs)
    debugSQL
}


### Use dbQuery(conn, SQL, 1) instead of dbQuery(conn, SQL)[[1]],
### it's much safer!
dbEasyQuery <- function(conn, SQL, j0=NA)
{
    if (get("debugSQL", envir=RTobjs)) {
        if (!is.character(SQL) || length(SQL) != 1L || is.na(SQL))
            stop("[debugSQL] 'SQL' must be a single string")
        cat("[debugSQL] SQL query: ", SQL, "\n", sep="")
        st <- system.time(data0 <- dbGetQuery(conn, SQL))
        cat("[debugSQL]      time: ", st["user.self"], " seconds\n", sep="")
    } else {
        data0 <- dbGetQuery(conn, SQL)
    }
    if (is.na(j0))
        return(data0)
    ## Needed to deal properly with data frame with 0 column ("NULL data
    ## frames with 0 rows") returned by RSQLite when the result of a SELECT
    ## query has 0 row
    if (nrow(data0) == 0L)
        character(0)
    else
        data0[[j0]]
}

### TODO: Put this in AnnotationDbi.
queryAnnotationDb <- function(annotationdb, sql)
{
    AnnotationDbi:::dbEasyQuery(dbconn(annotationdb),
                                paste(sql, collapse="\n"))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Data frame related.
###
### TODO: Find a better home for these low-level data.frame utils.
###

### Not data.frame specific. Would also work on any matrix-like object.
has_col <- function(x, colnames) {colnames %in% colnames(x)}

makeZeroRowDataFrame <- function(col2class)
{
    if (!is.character(col2class) || is.null(names(col2class)))
        stop("'col2class' must be a named character vector")
    as.data.frame(lapply(col2class, function(class) get(class)()),
                  stringsAsFactors=FALSE)
}

### Sets the class of (all or some of) the columns of a data.frame.
### Typical use:
###   x <- setDataFrameColClass(x, c(colA="integer", colB="factor"))
### Note that if 'x' has more than one "colA" col, then *all* of them are
### coerced to integer.
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
    if (drop.extra.cols) {
        col_idx <- which(colnames(x) %in% names(col2class))
    } else {
        col_idx <- seq_len(ncol(x))
    }
    tmp <- lapply(col_idx,
                  function(j)
                  {
                      col <- x[[j]]
                      colname <- colnames(x)[j]
                      if (!(colname %in% names(col2class)))
                          return(col)
                      class <- col2class[[colname]]
                      FUNname <- paste("as", class, sep=".")
                      if (exists(FUNname) && is.function(FUN <- get(FUNname)))
                          return(FUN(col))
                      as(col, class)
                  })
    names(tmp) <- colnames(x)[col_idx]
    return(data.frame(tmp, check.names=FALSE, stringsAsFactors=FALSE))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ID assignment and/or reassignment.
###

### Returns the vector of ids such that 'unique(x)[ids, ]' is identical
### to 'x' (in the same way that 'levels(f)[f]' is identical to
### 'as.vector(f)' when 'f' is a character factor).
### This unambiguously defines 'ids'. In particular, it's not Locale
### specific, despite the fact that the current implementation uses a
### sorting approach.
### TODO: Remove! (not used anymore)
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

