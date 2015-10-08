### =========================================================================
### The transcripts(), exons(), cds() and promoters() extractors
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers for building SQL queries
###

.build_SQL_FROM <- function(joins, join_type="INNER")
{
    joins_len <- length(joins)
    stopifnot(joins_len %% 2L == 1L)
    SQL <- joins[[1L]]
    if (joins_len != 1L) {
        ON_idx <- 2L * seq_len(joins_len %/% 2L)
        ON <- joins[ON_idx]
        Rtables <- joins[ON_idx + 1L]
        SQL <- c(SQL, paste0(join_type, " JOIN ", Rtables, " ON (", ON, ")"))
    }
    SQL
}

.build_SQL_WHERE <- function(vals)
{
    if (length(vals) == 0L)
        return("")
    sql <-
      lapply(seq_len(length(vals)), function(i) {
               v <- vals[[i]]
               if (!is.numeric(v))
                 v <- paste0("'", v, "'")
               v <- paste0("(", paste0(v, collapse=","), ")")
               v <- paste0(names(vals)[i], " IN ", v)
               paste0("(", v, ")")
            })
    paste0(unlist(sql), collapse=" AND ")
}

.build_SQL_SELECT <- function(columns, joins, distinct=FALSE,
                              vals=list(), orderby=character(0))
{
    SQL <- "SELECT"
    if (distinct)
        SQL <- c(SQL, "DISTINCT")
    SQL <- c(SQL, paste0(columns, collapse=", "),
             "FROM", .build_SQL_FROM(joins))
    if (length(vals) != 0L)
        SQL <- c(SQL, "WHERE", .build_SQL_WHERE(vals))
    if (length(orderby) != 0L)
        SQL <- c(SQL, "ORDER BY", paste0(orderby, collapse=", "))
    SQL
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### DB schema
###

CHROMNFO_DESC <- list(
    cols=c("_chrom_id",
           "chrom",
           "length",
           "is_circular"),
    Pcol="_chrom_id"
)

TRANSCRIPT_DESC <- list(
    cols=c(id="_tx_id",
           name="tx_name",
           type="tx_type",
           chrom="tx_chrom",
           strand="tx_strand",
           start="tx_start",
           end="tx_end"),
    Pcol="_tx_id"
)

EXON_DESC <- list(
    cols=c(id="_exon_id",
           name="exon_name",
           chrom="exon_chrom",
           strand="exon_strand",
           start="exon_start",
           end="exon_end"),
    Pcol="_exon_id"
)

CDS_DESC <- list(
    cols=c(id="_cds_id",
           name="cds_name",
           chrom="cds_chrom",
           strand="cds_strand",
           start="cds_start",
           end="cds_end"),
    Pcol="_cds_id"
)

SPLICING_DESC <- list(
    cols=c("_tx_id", "exon_rank", "_exon_id", "_cds_id")
)

GENE_DESC <- list(
    cols=c("gene_id", "_tx_id")
)

### Order of tables matters! "transcript" must be before "splicing" and "gene",
### and "exon" and "cds" must be before "splicing". See .column2table() below
### why.
DB_DESC <- list(
    chrominfo=CHROMNFO_DESC,
    transcript=TRANSCRIPT_DESC,
    exon=EXON_DESC,
    cds=CDS_DESC,
    splicing=SPLICING_DESC,
    gene=GENE_DESC
)

### Tables "transcript", "exon", and "cds", must have these tags (at a minimum).
CORE_TAGS <- c("id", "chrom", "strand", "start", "end")

### The "splicing right tables" can be bundled to the "splicing" table with
### a LEFT JOIN using the SPLICING_JOIN_USING columns.
SPLICING_RTABLES <- c("transcript", "exon", "cds")
SPLICING_JOIN_USING <- setNames(c("_tx_id", "_exon_id", "_cds_id"),
                                SPLICING_RTABLES)

### When the same column belongs to more than one table (e.g. "_tx_id",
### "_exon_id", or "_cds_id"), then the table for which the column is a
### primary key is chosen.
.column2table <- function(columns)
{
    if (length(columns) == 0L)
        return(character(0))
    sapply(columns, function(column) {
        for (table in names(DB_DESC)) {
            if (column %in% DB_DESC[[table]]$cols)
                return(table)
        }
        stop(column, ": unknown db column")
    })
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .select_features()
###

.as_qualified <- function(tables, columns) paste(tables, columns, sep=".")

.join_splicing_Rtables <- function(tables=character(0))
{
    tables <- unique(tables)
    join_order <- c("splicing", SPLICING_RTABLES)
    if (!all(tables %in% join_order))
        stop("invalid tables")
    tables <- c("splicing", tables)
    ## Order tables & remove duplicates.
    tables <- intersect(join_order, tables)
    if (length(tables) == 1L)
        return(tables)
    joins <- character(2L * length(tables) - 1L)
    ON_idx <- 2L * seq_len(length(tables) - 1L)
    Rtables <- tables[-1L]
    USING <- SPLICING_JOIN_USING[Rtables]
    Lcolumns <- .as_qualified("splicing", USING)
    Rcolumns <- .as_qualified(Rtables, USING)
    ON <- paste(Lcolumns, Rcolumns, sep="=")
    joins[ON_idx] <- ON     
    joins[c(1L, ON_idx + 1L)] <- tables
    joins
}

.join_tables <- function(tables)
{
    tables <- unique(tables)
    if (length(tables) == 1L)
        return(tables)
    if (any(tables %in% c("exon", "cds")))
        tables <- c(tables, "splicing")
    ## Order tables & remove duplicates.
    join_order <- c("transcript", "gene", "splicing", "exon", "cds")
    tables <- intersect(join_order, tables)
    joins <- character(2L * length(tables) - 1L)
    ON_idx <- 2L * seq_len(length(tables) - 1L)
    ON <- sapply(2:length(tables), function(i) {
        Rtable <- tables[[i]]
        if (Rtable == "exon") {
            USING <- "_exon_id"
            Ltable <- "splicing"
        } else if (Rtable == "cds") {
            USING <- "_cds_id"
            Ltable <- "splicing"
        } else {
            USING <- "_tx_id"
            Ltable <- tables[[i - 1L]]
        }
        Lcolumn <- .as_qualified(Ltable, USING)
        Rcolumn <- .as_qualified(Rtable, USING)
        paste(Lcolumn, Rcolumn, sep="=")
    })
    joins[ON_idx] <- ON
    joins[c(1L, ON_idx + 1L)] <- tables
    joins
}

.select_features <- function(txdb, columns=character(0), vals=list(),
                             ptable, core_columns)
{
    pkey <- orderby <- core_columns[["id"]]  # primary key
    names(columns) <- .column2table(columns)

    ## 1st SQL query: extract stuff from the primary table.
    columns1 <- union(core_columns, columns[names(columns) == ptable])
    where_columns <- names(vals)
    where_tables <- .column2table(where_columns)
    joins <- .join_tables(c(ptable, where_tables))
    if (length(joins) != 1L) {
        columns1 <- .as_qualified(.column2table(columns1), columns1)
        names(vals) <- .as_qualified(where_tables, where_columns)
        orderby <- .as_qualified(.column2table(orderby), orderby)
    }
    SQL <- .build_SQL_SELECT(columns1, joins, distinct=TRUE,
                             vals=vals, orderby=orderby)
    df1 <- queryAnnotationDb(txdb, SQL)

    ## Additional SQL queries: 1 additional query per secondary table with the
    ## exception that "splicing right tables" are treated as the virtual single
    ## table obtained by LEFT JOIN'ing them together.
    foreign_columns <- columns[names(columns) != ptable]
    splicing_bundle <- c("splicing", SPLICING_RTABLES)
    bundle_idx <- names(foreign_columns) %in% splicing_bundle
    names(foreign_columns)[bundle_idx] <- "splicing"
    foreign_columns <- split(foreign_columns, names(foreign_columns))
    secondary_tables <- names(foreign_columns)
    names(secondary_tables) <- secondary_tables
    df_list <- lapply(secondary_tables, function(table) {
        columns2 <- foreign_columns[[table]]
        if (length(vals) == 0L) {
            vals2 <- list()
        } else {
            vals2 <- setNames(list(df1[[pkey]]), pkey)
        }
        if (table == "splicing") {
            tables2 <- .column2table(columns2)
            joins <- .join_splicing_Rtables(tables2)
            from <- paste0(.build_SQL_FROM(joins, "LEFT"), collapse=" ")
            columns2 <- .as_qualified(tables2, columns2)
            columns2 <- c(pkey, columns2)
            if (ptable == "transcript") {
                orderby <- c(pkey, "exon_rank")
            } else {
                orderby <- c(pkey, "splicing._tx_id")
            }
            SQL <- .build_SQL_SELECT(columns2, from,
                                     vals=vals2, orderby=orderby)
        } else {
            joins <- .join_tables(c(ptable, table))
            columns2 <- c(pkey, columns2)
            columns2 <- .as_qualified(.column2table(columns2), columns2)
            names(vals2) <- .as_qualified(.column2table(names(vals2)),
                                          names(vals2))
            SQL <- .build_SQL_SELECT(columns2, joins, distinct=TRUE,
                                     vals=vals2)
        }
        queryAnnotationDb(txdb, SQL)
    })
    c(setNames(list(df1), ptable), df_list)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_features_as_GRanges()
###

.make_DataFrame_from_df_list <- function(df_list)
{
    DF1 <- DataFrame(df_list[[1L]], check.names=FALSE)
    if (length(df_list) == 1L)
        return(DF1)
    pkey <- names(DF1)[[1L]]
    ids <- DF1[[pkey]]
    DF_list <- lapply(2:length(df_list), function(i) {
        df <- df_list[[i]]
        stopifnot(identical(names(df)[[1L]], pkey))
        f <- factor(df[[pkey]], levels=ids)
        DataFrame(lapply(df[-1L], function(col) unname(splitAsList(col, f))),
                  check.names=FALSE)
    })
    do.call(DataFrame, c(list(DF1), DF_list, list(check.names=FALSE)))
}

.as_db_columns <- function(columns)
    sub("^(tx_id|exon_id|cds_id)$", "_\\1", columns)

.extract_features_as_GRanges <- function(txdb, ptable,
                                         columns=character(0), vals=list())
{
    db_columns <- .as_db_columns(columns)
    names(vals) <- .as_db_columns(names(vals))
    core_columns <- DB_DESC[[ptable]]$cols[CORE_TAGS]
    df_list <- .select_features(txdb, db_columns, vals, ptable, core_columns)
    DF <- .make_DataFrame_from_df_list(df_list)
    ans <- makeGRangesFromDataFrame(DF, seqinfo=get_TxDb_seqinfo0(txdb))
    mcols(ans) <- setNames(DF[db_columns], columns)
    keep_user_seqlevels_from_TxDb(ans, txdb)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Primary extractors: transcripts(), exons(), cds(), and genes().
###
### TODO: Rename the 'vals' arg -> 'filter' so it's consistent with the
### 'filters' arg of makeTxDbFromBiomart() which we should also rename
### 'filter'. Also place it *after* the 'columns' argument.
### Other proposal: rename the 'columns' arg -> 'colnames'
###

.dbSchemaHasTxType <- function(txdb){
    txFields <- dbListFields(dbconn(txdb),"transcript")
    "tx_type" %in% txFields
}

## For converting user arguments FROM the UC style down to what we use
## internally
translateCols <- function(columns, txdb){
    ## preserve any names
    oriColNames <- names(columns)
    ## and save the original column strings
    oriCols <- columns
    
    oriLen <- length(columns) ## not always the same as length(oriCols)
    ## get the available abbreviations as a translation vector (exp)
    names <- .makeColAbbreviations(txdb)
    if(("TXTYPE" %in% columns) && !.dbSchemaHasTxType(txdb)){
        names <- c(names, c(tx_type='TXTYPE'))
    }
    exp <- sub("^_","", names(names))
    names(exp) <- names

    ## Then replace only those IDs that match the UC names
    m <- match(oriCols, names(exp))
    idx = which(!is.na(m))
    columns[idx] <- exp[m[idx]]
    
    if(length(columns) == oriLen && is.null(oriColNames)){
        names(columns) <- oriCols 
    }else if(length(columns) == oriLen && !is.null(oriColNames)){
        names(columns) <- oriColNames         
    }else if(length(columns) != oriLen){
        stop("names were lost in translateCols() helper")
    }    
    columns
}

## me a named list from the metadata data.frame
.makeMetadataList <- function(meta){
    lst <- as.list(meta[,2])
    names(lst) <- meta[,1]
    lst
}

## assign this to the metadata list in relevant object
.assignMetadataList <- function(obj, txdb){
    metadata(obj)[[1]] <- .makeMetadataList(metadata(txdb))
    names(metadata(obj))[[1]] <- "genomeInfo"
    obj
}

.extractFromTxDb <- function(txdb, ptable, columns=character(0), vals=NULL)
{
    user_columns <- columns
    columns <- translateCols(columns, txdb)
    if (is.null(vals))
        vals <- list()
    names(vals) <- translateCols(names(vals), txdb)
    ans <- .extract_features_as_GRanges(txdb, ptable, columns, vals)
    names(mcols(ans)) <- if (is.null(names(user_columns))) user_columns
                         else names(user_columns)
    .assignMetadataList(ans, txdb)
}

setGeneric("transcripts", function(x, ...) standardGeneric("transcripts"))

setMethod("transcripts", "TxDb",
    function(x, vals=NULL, columns=c("tx_id", "tx_name"))
        .extractFromTxDb(x, "transcript", columns=columns, vals=vals)
)

setGeneric("exons", function(x, ...) standardGeneric("exons"))

setMethod("exons", "TxDb",
    function(x, vals=NULL, columns="exon_id")
        .extractFromTxDb(x, "exon", columns=columns, vals=vals)
)

setGeneric("cds", function(x, ...) standardGeneric("cds"))

setMethod("cds", "TxDb",
    function(x, vals=NULL, columns="cds_id")
        .extractFromTxDb(x, "cds", columns=columns, vals=vals)
)

setGeneric("genes", function(x, ...) standardGeneric("genes"))

.relist_col <- function(x, skeleton)
{
   if (is.list(x) || (is(x, "List") && !is(x, "Ranges")))
       return(IRanges:::regroupBySupergroup(x, skeleton))
   relist(x, skeleton)
}

### Used in SplicingGraphs! Move this to IRanges (and check similarity with
### "collapse" method for DataFrame objects).
.collapse_df <- function(df, skeleton)
{
    ## FIXME: endoapply() on a DataFrame object is broken when applying
    ## a function 'FUN' that modifies the nb of rows. Furthermore, the
    ## returned object passes validation despite being broken! Fix it
    ## in IRanges.
    ans <- endoapply(df, function(x) unique(.relist_col(x, skeleton)))
    ## Fix the broken DataFrame returned by endoapply().
    ans@nrows <- length(skeleton)
    ans@rownames <- NULL
    ans
}

### If 'single.strand.genes.only' is TRUE (the default), then genes that
### have exons located on both strands of the same chromosome, or on 2
### different chromosomes are dropped. In that case, the genes are returned
### in a GRanges object. Otherwise, they're returned in a GRangesList object
### with the metadata columns requested thru 'columns' set at the top level.
.TxDb.genes <- function(x, vals=NULL, columns="gene_id",
                        single.strand.genes.only=TRUE)
{
    if (!is.character(columns))
        stop("'columns' must be a character vector")
    if (!isTRUEorFALSE(single.strand.genes.only))
        stop("'single.strand.genes.only' must be TRUE or FALSE")
    columns2 <- union(columns, "gene_id")
    tx <- transcripts(x, vals=vals, columns=columns2)

    ## Unroll 'tx' along the 'gene_id' metadata column.
    ## Note that the number of genes per transcript will generally be 1 or 0.
    ## But we also want to handle the situation where it's > 1 which happens
    ## when the same transcript is linked to more than 1 gene (because this
    ## may happen one day and is the reason behind the choice to represent
    ## the 'gene_id' as a CharacterList object instead of a character vector).
    gene_id <- mcols(tx)[ , "gene_id"]
    ngene_per_tx <- elementLengths(gene_id)
    tx <- tx[rep.int(seq_along(ngene_per_tx), ngene_per_tx)]
    mcols(tx)$gene_id <- unlist(gene_id, use.names=FALSE)

    ## Split 'tx' by gene.
    tx_by_gene <- split(tx, mcols(tx)$gene_id)

    ## Turn inner mcols into outter mcols by relist()'ing them.
    inner_mcols <- mcols(tx_by_gene@unlistData)[columns]
    mcols(tx_by_gene@unlistData) <- NULL
    outter_mcols <- .collapse_df(inner_mcols, tx_by_gene)
    gene_id <- outter_mcols$gene_id
    if (!is.null(gene_id)) {
        stopifnot(identical(names(tx_by_gene), as.character(gene_id)))
        outter_mcols$gene_id <- names(tx_by_gene)
    }
    mcols(tx_by_gene) <- outter_mcols

    ## Compute the gene ranges.
    genes <- range(tx_by_gene)

    if (!single.strand.genes.only)
        return(genes)

    keep_idx <- which(elementLengths(genes) == 1L)
    genes <- genes[keep_idx]
    ans <- unlist(genes, use.names=FALSE)
    mcols(ans) <- mcols(genes)
    names(ans) <- names(genes)
    ans
}

setMethod("genes", "TxDb", .TxDb.genes)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "promoters" method
###
### generic is in IRanges
###

setMethod("promoters", "TxDb",
    function(x, upstream=2000, downstream=200, ...)
    {
        gr <- transcripts(x, ...)
        promoters(gr, upstream=upstream, downstream=downstream)
    }
) 

