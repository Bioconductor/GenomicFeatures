### =========================================================================
### Miscellaneous low-level utils
### -------------------------------------------------------------------------
###
### Unless stated otherwise, nothing in this file is exported.
###


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

dbEasyPreparedQuery <- function(conn, SQL, bind.data)
{
    ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
    ## when the nb of rows to insert is 0, hence the early bail out.
    if (nrow(bind.data) == 0L)
        return()
    if (get("debugSQL", envir=RTobjs)) {
        if (!is.character(SQL) || length(SQL) != 1L || is.na(SQL))
            stop("[debugSQL] 'SQL' must be a single string")
        cat("[debugSQL] SQL prepared query: ", SQL, "\n", sep="")
        cat("[debugSQL]     dim(bind.data): ",
            paste(dim(bind.data), collapse=" x "), "\n", sep="")
        st <- system.time({
                  dbBegin(conn)
                  dbGetQuery(conn, SQL, unclass(unname(bind.data)))
                  dbCommit(conn)})
        cat("[debugSQL]               time: ", st["user.self"],
            " seconds\n", sep="")
    } else {
        dbBegin(conn)
        dbGetQuery(conn, SQL, unclass(unname(bind.data)))
        dbCommit(conn)
    }
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

check_colnames <- function(x, required_colnames,
                              optional_colnames=character(0),
                              x_what="x")
{
    if (!is.data.frame(x))
        stop(wmsg("'", x_what, "' must be a data frame"))
    x_colnames <- colnames(x)
    if (!all(required_colnames %in% x_colnames))
        stop(wmsg("'", x_what, "' must have at least the following col",
                  ifelse(length(required_colnames) >= 2L, "s", ""),
                  ": ", paste(required_colnames, collapse=", ")))
    supported_colnames <- c(required_colnames, optional_colnames)
    is_supported_col <- x_colnames %in% supported_colnames
    if (any(duplicated(x_colnames[is_supported_col])))
        stop(wmsg("'", x_what, "' has duplicated colnames"))
    ignored_idx <- which(!is_supported_col)
    if (length(ignored_idx) != 0L)
        warning(wmsg("ignoring the following col",
                     ifelse(length(ignored_idx) >= 2L, "s", ""),
                     " in '", x_what, "': ",
                     paste(x_colnames[ignored_idx], collapse=", ")))
}

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
    if (!isSingleString(join_colname) || !has_col(x, join_colname))
        stop("'join_colname' must be a valid colname for 'x'")
    if (!is.vector(name2val) && !is.factor(name2val))
        stop("'name2val' must be a vector (or factor)")
    if (!is.atomic(name2val) || is.null(names(name2val)))
        stop("'name2val' must be atomic and have names")
    if (!isSingleString(vals_colname))
        stop("invalid 'vals_colname'")
    join_col <- as.character(x[[join_colname]])
    common_names <- intersect(join_col, names(name2val))
    name2val <- name2val[names(name2val) %in% common_names]
    x <- S4Vectors:::extract_data_frame_rows(x, join_col %in% common_names)
    tmp <- split(as.vector(name2val), names(name2val))
    ## as.character() is required below just because 'x[[join_colname]]'
    ## could be a factor (subsetting by a factor is equivalent to subsetting
    ## by an integer vector but this is not what we want here).
    tmp <- tmp[as.character(x[[join_colname]])]
    x <- S4Vectors:::extract_data_frame_rows(x,
                             rep.int(seq_len(nrow(x)), elementNROWS(tmp)))
    if (nrow(x) == 0L)
        vals <- name2val[FALSE]
    else if (is.factor(name2val))
        vals <- factor(unname(unlist(tmp)), levels=levels(name2val))
    else
        vals <- unname(unlist(tmp))
    x[[vals_colname]] <- vals
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ID assignment and/or reassignment.
###

### Mimicking the interface of chartr().
### If 'old_ids' and 'ids' are character vectors, then
### 'translateIds(old_ids, new_ids, ids)' is equivalent to
### 'names(new_ids) <- old_ids; new_ids[ids]'.
translateIds <- function(old_ids, new_ids, ids)
{
    if (!is.atomic(old_ids) || !is.atomic(new_ids) || !is.atomic(ids))
        stop("'old_ids', 'new_ids' and 'ids' must be atomic vectors")
    if (length(old_ids) != length(new_ids))
        stop("'old_ids' and 'new_ids' must have the same length")
    new_ids[match(ids, old_ids)]
}

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

### 'name' and 'type' each must be either NULL, or an integer or character
### vector, or a factor. Possibly with NAs. One of them can be NULL but not
### both of them. If none of them is NULL, they must have the same length N.
### Returns an integer vector of length N.
.rank_name_type <- function(name, type)
{
    ## We set LC_COLLATE to C so (1) the ouput of rank() on a character vector
    ## is platform/country independent, and (2) it will behave the same way
    ## when called in the context of the unit tests run by 'R CMD check'
    ## ('R CMD check' also sets the LC_COLLATE to C when running the tests)
    ## vs when called in the context of an interactive session.
    ## TODO: Maybe we should have a strrank() function in
    ## S4Vectors/R/str-utils.R for portable/deterministic ranking of a
    ## character vector.
    prev_locale <- Sys.getlocale("LC_COLLATE")
    Sys.setlocale("LC_COLLATE", "C")
    on.exit(Sys.setlocale("LC_COLLATE", prev_locale))
    if (!is.integer(name)) {
        if (is.null(name)) {
            name <- integer(length(type))
        } else if (is.numeric(name)) {
            name <- as.integer(name)
        } else {
            if (!is.character(name))
                name <- as.character(name)
            name <- rank(name, na.last="keep", ties.method="min")
        }
    }
    ## Features with no name (e.g. tx_name is NA) go last.
    name[is.na(name)] <- .Machine$integer.max
    if (!is.integer(type)) {
        if (is.null(type)) {
            type <- integer(length(name))
        } else if (is.factor(type) || is.numeric(type)) {
            type <- as.integer(type)
        } else {
            if (!is.character(type))
                type <- as.character(type)
            type <- rank(type, na.last="keep", ties.method="min")
        }
    }
    ## Features with no type (e.g. tx_type is NA) go last.
    type[is.na(type)] <- .Machine$integer.max
    oo <- orderIntegerPairs(name, type)
    ans <- integer(length(oo))
    ans[oo] <- seq_along(oo)
    sm <- selfmatchIntegerPairs(name, type)
    ans[sm]
}

### 'chrom_ids' (integer vector, no NAs), 'strand' (character vector, factor,
### or anything supported by a "strand" method, no NAs), 'start' (integer
### vector, no NAs), and 'end' (integer vector, no NAs) must have the same
### length N (number of features).
### If 'name' is not NULL, it must be character vector or factor of length N,
### possibly with NAs.
### Returns an integer vector of length N containing one id per feature.
makeFeatureIds <- function(chrom_ids, strand, start, end,
                           name=NULL, type=NULL,
                           same.id.for.dups=FALSE)
{
    if (is.factor(strand)) {
        ## If levels contain "+", "-" and/or "*" in the wrong order then
        ## we coerce back to character.
        m <- match(levels(strand), levels(strand()))
        m <- m[!is.na(m)]
        if (!all(diff(m) >= 1L))
            strand <- as.character(strand)
    }
    if (!is.factor(strand))
        strand <- strand(strand)
    a <- chrom_ids
    b <- as.integer(strand)
    c <- start
    d <- end
    if (!(is.null(name) && is.null(type))) {
        a <- 3L * a + b
        b <- c
        c <- d
        d <- .rank_name_type(name, type)
    }
    if (!same.id.for.dups) {
        oo <- orderIntegerQuads(a, b, c, d)
        ans <- integer(length(oo))
        ans[oo] <- seq_len(length(oo))
        return(ans)
    }
    ## There should be a better way to do this...
    is_not_dup <- !duplicatedIntegerQuads(a, b, c ,d)
    ua <- a[is_not_dup]
    ub <- b[is_not_dup]
    uc <- c[is_not_dup]
    ud <- d[is_not_dup]
    oo <- orderIntegerQuads(ua, ub, uc, ud)
    ua <- ua[oo]
    ub <- ub[oo]
    uc <- uc[oo]
    ud <- ud[oo]
    matchIntegerQuads(a, b, c, d, ua, ub, uc, ud)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Miscellaneous.
###

### Global character vector to hold default names for circular sequences.
### This is exported!
DEFAULT_CIRC_SEQS <- c(
    ## Mitochondrial genome
    "chrM", "MT", "MtDNA", "mit", "Mito", "mitochondrion",
    "dmel_mitochondrion_genome",
    ## Chloroplast genome
    "Pltd", "ChrC", "Pt", "chloroplast", "Chloro",
    ## Plasmid (yeast)
    "2micron", "2-micron", "2uM"
)

### AFAIK UCSC doesn't flag circular sequences.
### As of Sep 21, 2010 (Ensembl release 59), Ensembl was still not flagging
### circular sequences in their db (see this thread for the details
### http://lists.ensembl.org/pipermail/dev/2010-September/000139.html),
make_circ_flags_from_circ_seqs <- function(seqlevels,
                                           circ_seqs=DEFAULT_CIRC_SEQS)
{
    if (!is.character(seqlevels))
        stop(wmsg("'seqlevels' must be a character vector"))
    if (identical(circ_seqs, DEFAULT_CIRC_SEQS)) {
        ## 'circ_seqs' is set to the default (very likely the user did NOT
        ## specify this argument).
        seqlevels <- tolower(seqlevels)
        circ_seqs <- tolower(circ_seqs)
        circ_flags <- rep.int(NA, length(seqlevels))
        circ_flags[seqlevels %in% circ_seqs] <- TRUE
    } else {
        ## The user specified the 'circ_seqs' argument.
        if (!is.character(circ_seqs) ||
            any(circ_seqs %in% c(NA_character_, "")))
            stop(wmsg("'circ_seqs' must be a character vector with no NAs ",
                      "and no empty strings"))
        bad_circ_seqs <- setdiff(circ_seqs, seqlevels)
        if (length(bad_circ_seqs) != 0L) {
            in1string <- paste0(bad_circ_seqs, collapse=", ")
            stop(wmsg("The following chromosome names in 'circ_seqs' are ",
                      "not found in the TxDb object to be made: ", in1string))
        }
        circ_flags <- seqlevels %in% circ_seqs
    }
    circ_flags
}

### 'exon_count' must be a vector of positive integers and 'tx_strand' a
### character vector with "+" or "-" values. Both vectors must have the
### same length.
makeExonRankCol <- function(exon_count, tx_strand)
{
    ans <- lapply(seq_len(length(exon_count)),
        function(i)
        {
            if (tx_strand[i] == "+")
                seq_len(exon_count[i])
            else
                (exon_count[i]):1L
        }
    )
    unlist(ans)
}
