### =========================================================================
### The transcripts(), exons(), cds() and promoters() extractors
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers for building SQL queries
###

.tables_in_joins <- function(joins)
{
    joins_len <- length(joins)
    stopifnot(joins_len %% 2L == 1L)
    joins[seq(1L, joins_len, by=2L)]
}

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
SPLICING_BUNDLE <- c("splicing", SPLICING_RTABLES)

### When the same column belongs to more than one table (e.g. "_tx_id",
### "_exon_id", or "_cds_id"), then the table for which the column is a
### primary key is chosen by default. This behavior can be changed by passing
### the name of a table to 'from_table' in which case the priority is given to
### that table.
.column2table <- function(columns, from_table=NA)
{
    if (length(columns) == 0L)
        return(character(0))
    tables <- sapply(columns,
        function(column) {
            for (table in names(DB_DESC)) {
                if (column %in% DB_DESC[[table]]$cols)
                    return(table)
            }
            stop(column, ": unknown db column")
        }
    )
    if (!is.na(from_table))
        tables[columns %in% DB_DESC[[from_table]]$cols] <- from_table
    tables
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .select_from_INNER_JOIN() and .select_from_SPLICING_BUNDLE()
###

.as_qualified <- function(tables, columns) paste(tables, columns, sep=".")

.join_tables <- function(tables)
{
    tables <- unique(tables)
    if (length(tables) == 1L)
        return(tables)
    if (any(tables %in% c("exon", "cds")))
        tables <- c(tables, "splicing")
    ## Order tables & remove duplicates.
    join_order <- c("transcript", "splicing", "exon", "cds", "gene")
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
            Ltable <- tables[[1L]]
        }
        Lcolumn <- .as_qualified(Ltable, USING)
        Rcolumn <- .as_qualified(Rtable, USING)
        paste(Lcolumn, Rcolumn, sep="=")
    })
    joins[ON_idx] <- ON
    joins[c(1L, ON_idx + 1L)] <- tables
    joins
}

.join_splicing_Rtables <- function(tables=character(0))
{
    if (!all(tables %in% SPLICING_BUNDLE))
        stop("all tables must be in SPLICING_BUNDLE")
    tables <- c("splicing", tables)
    ## Order tables & remove duplicates.
    tables <- intersect(SPLICING_BUNDLE, tables)
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

### The columns in 'columns' + those involved thru 'vals' and 'orderby' are
### collected and their corresponding tables are INNER JOIN'ed.
.select_from_INNER_JOIN <- function(txdb, table, columns, vals=list(),
                                    orderby=character(0))
{
    tables <- .column2table(columns, from_table=table)
    where_columns <- names(vals)
    where_tables <- .column2table(where_columns, from_table=table)
    joins <- .join_tables(c(table, tables, where_tables))
    orderby_tables <- .column2table(orderby, from_table=table)
    stopifnot(all(orderby_tables %in% .tables_in_joins(joins)))
    use_joins <- length(joins) > 1L
    if (use_joins) {
        columns <- .as_qualified(tables, columns)
        names(vals) <- .as_qualified(where_tables, where_columns)
        orderby <- .as_qualified(orderby_tables, orderby)
    }
    ## .build_SQL_SELECT() generates an INNER JOIN.
    SQL <- .build_SQL_SELECT(columns, joins, distinct=use_joins,
                             vals=vals, orderby=orderby)
    queryAnnotationDb(txdb, SQL)
}

### Can only involve columns (thru 'columns', 'vals', and 'orderby') that
### belong to the tables in SPLICING_BUNDLE at the moment.
.select_from_SPLICING_BUNDLE <- function(txdb, columns,
                                         vals=list(), orderby=character(0))
{
    tables <- .column2table(columns, from_table="splicing")
    where_columns <- names(vals)
    where_tables <- .column2table(where_columns, from_table="splicing")
    orderby_tables <- .column2table(orderby, from_table="splicing")
    joins <- .join_splicing_Rtables(c(tables, where_tables, orderby_tables))
    use_joins <- length(joins) > 1L
    if (use_joins) {
        columns <- .as_qualified(tables, columns)
        names(vals) <- .as_qualified(where_tables, where_columns)
        orderby <- .as_qualified(orderby_tables, orderby)
    }
    ## .build_SQL_SELECT() would generate an INNER JOIN but we want a
    ## LEFT JOIN.
    from <- paste0(.build_SQL_FROM(joins, "LEFT"), collapse=" ")
    SQL <- .build_SQL_SELECT(columns, from, distinct=FALSE,
                             vals=vals, orderby=orderby)
    queryAnnotationDb(txdb, SQL)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_features()
###

.extract_features <- function(txdb, proxy_table, columns=character(0),
                              vals=list(), core_columns)
{
    proxy_column <- orderby <- core_columns[["id"]]
    names(columns) <- .column2table(columns, from_table=proxy_table)

    ## 1st SQL query: extract stuff from the proxy table.
    columns1 <- union(core_columns, columns[names(columns) == proxy_table])
    df1 <- .select_from_INNER_JOIN(txdb, proxy_table, columns1, vals=vals,
                                   orderby=orderby)

    ## Additional SQL queries: 1 additional query per satellite table with the
    ## exception that the satellite tables that belong to SPLICING_BUNDLE are
    ## treated as the virtual single table obtained by LEFT JOIN'ing them
    ## together.
    foreign_columns <- columns[names(columns) != proxy_table]
    bundle_idx <- names(foreign_columns) %in% SPLICING_BUNDLE
    names(foreign_columns)[bundle_idx] <- "splicing"
    foreign_columns <- split(foreign_columns, names(foreign_columns))
    satellite_tables <- names(foreign_columns)
    names(satellite_tables) <- satellite_tables
    df_list <- lapply(satellite_tables, function(satellite_table) {
        columns2 <- foreign_columns[[satellite_table]]
        if (length(vals) == 0L) {
            vals2 <- list()
        } else {
            vals2 <- setNames(list(df1[[proxy_column]]), proxy_column)
        }
        if (satellite_table == "splicing") {
            columns2 <- c(proxy_column, columns2)
            if (proxy_table == "transcript") {
                orderby <- c(proxy_column, "exon_rank")
            } else {
                orderby <- c(proxy_column, "_tx_id")
            }
            .select_from_SPLICING_BUNDLE(txdb, columns2,
                                         vals=vals2, orderby=orderby)
        } else if (satellite_table == "gene") {
            columns2 <- c(proxy_column, columns2)
            orderby <- c("_tx_id", "gene_id")
            .select_from_INNER_JOIN(txdb, "gene", columns2,
                                    vals=vals2, orderby=orderby)
        } else {
            stop(satellite_table, ": unsupported satellite table")
        }
    })
    c(setNames(list(df1), proxy_table), df_list)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_features_as_GRanges()
###

.make_DataFrame_from_df_list <- function(df_list)
{
    DF1 <- DataFrame(df_list[[1L]], check.names=FALSE)
    if (length(df_list) == 1L)
        return(DF1)
    proxy_column <- names(DF1)[[1L]]
    ids <- DF1[[proxy_column]]
    DF_list <- lapply(2:length(df_list), function(i) {
        df <- df_list[[i]]
        stopifnot(identical(names(df)[[1L]], proxy_column))
        f <- factor(df[[proxy_column]], levels=ids)
        DataFrame(lapply(df[-1L], function(col) unname(splitAsList(col, f))),
                  check.names=FALSE)
    })
    do.call(DataFrame, c(list(DF1), DF_list, list(check.names=FALSE)))
}

.as_db_columns <- function(columns)
    sub("^(tx_id|exon_id|cds_id)$", "_\\1", columns)

.extract_features_as_GRanges <- function(txdb, proxy_table,
                                         columns=character(0), vals=list())
{
    db_columns <- .as_db_columns(columns)
    names(vals) <- .as_db_columns(names(vals))
    core_columns <- DB_DESC[[proxy_table]]$cols[CORE_TAGS]
    df_list <- .extract_features(txdb, proxy_table, db_columns,
                                 vals, core_columns)
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

.extractFromTxDb <- function(txdb, proxy_table, columns=character(0), vals=NULL)
{
    user_columns <- columns
    columns <- translateCols(columns, txdb)
    if (is.null(vals))
        vals <- list()
    names(vals) <- translateCols(names(vals), txdb)
    ans <- .extract_features_as_GRanges(txdb, proxy_table, columns, vals)
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

