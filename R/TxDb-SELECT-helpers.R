### =========================================================================
### Helpers for SELECT'ing stuff from a TxDb object
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level helpers (schema agnostic) for building SQL queries
###

.as_qualified <- function(tables, columns) paste(tables, columns, sep=".")

.tables_in_joins <- function(joins)
{
    joins_len <- length(joins)
    stopifnot(joins_len %% 2L == 1L)
    joins[seq(1L, joins_len, by=2L)]
}

### 'join_type' will be recycled to the nb of joins (= length(joins) %/% 2).
.build_SQL_FROM <- function(joins, join_type="INNER")
{
    joins_len <- length(joins)
    stopifnot(joins_len %% 2L == 1L)
    SQL <- joins[[1L]]
    if (joins_len == 1L)
        return(SQL)
    njoin <- joins_len %/% 2L
    stopifnot(length(join_type) == 1L || length(join_type) == njoin)
    ON_idx <- 2L * seq_len(njoin)
    ON <- joins[ON_idx]
    Rtables <- joins[ON_idx + 1L]
    c(SQL, paste0(join_type, " JOIN ", Rtables, " ON (", ON, ")"))
}

.build_SQL_FROM_splicing <- function(joins, cds_join_type="LEFT")
{
    joins_len <- length(joins)
    stopifnot(joins_len %% 2L == 1L)
    SQL <- joins[[1L]]
    if (joins_len == 1L)
        return(SQL)
    njoin <- joins_len %/% 2L
    join_type <- rep.int("INNER", njoin)
    if (joins[[length(joins)]] == "cds")
        join_type[[length(join_type)]] <- cds_join_type
    paste0(.build_SQL_FROM(joins, join_type), collapse=" ")
}

.build_SQL_WHERE <- function(filter)
{
    if (length(filter) == 0L)
        return("")
    sql <- lapply(seq_len(length(filter)),
             function(i) {
               fi <- filter[[i]]
               if (!is.numeric(fi))
                 fi <- paste0("'", fi, "'")
               fi <- paste0("(", paste0(fi, collapse=","), ")")
               fi <- paste0(names(filter)[i], " IN ", fi)
               paste0("(", fi, ")")
           })
    paste0(unlist(sql), collapse=" AND ")
}

.build_SQL_SELECT <- function(columns, joins, distinct=FALSE,
                              filter=list(), orderby=character(0))
{
    SQL <- "SELECT"
    if (distinct)
        SQL <- c(SQL, "DISTINCT")
    SQL <- c(SQL, paste0(columns, collapse=", "),
             "FROM", .build_SQL_FROM(joins))
    if (length(filter) != 0L)
        SQL <- c(SQL, "WHERE", .build_SQL_WHERE(filter))
    if (length(orderby) != 0L)
        SQL <- c(SQL, "ORDER BY", paste0(orderby, collapse=", "))
    SQL
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .TXDB_join_tables() and .TXDB_join_splicing_Rtables()
###

.TXDB_join_tables <- function(tables)
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

.TXDB_join_splicing_Rtables <- function(tables=character(0))
{
    if (!all(tables %in% TXDB_SPLICING_BUNDLE))
        stop("all tables must be in TXDB_SPLICING_BUNDLE")
    tables <- c("splicing", tables)
    ## Order tables & remove duplicates.
    tables <- intersect(TXDB_SPLICING_BUNDLE, tables)
    if (length(tables) == 1L)
        return(tables)
    joins <- character(2L * length(tables) - 1L)
    ON_idx <- 2L * seq_len(length(tables) - 1L)
    Rtables <- tables[-1L]
    USING <- TXDB_SPLICING_JOIN_USING[Rtables]
    Lcolumns <- .as_qualified("splicing", USING)
    Rcolumns <- .as_qualified(Rtables, USING)
    ON <- paste(Lcolumns, Rcolumns, sep="=")
    joins[ON_idx] <- ON
    joins[c(1L, ON_idx + 1L)] <- tables
    joins
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TxDb_schema_version()
###

TxDb_schema_version <- function(txdb)
{
    version <- AnnotationDbi:::.getMetaValue(dbconn(txdb), "DBSCHEMAVERSION")
    numeric_version(version)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 2 flexible helpers for SELECT'ing stuff from a TxDb object:
###   - TxDb_SELECT_from_INNER_JOIN()
###   - TxDb_SELECT_from_splicing_bundle()
### They should satisfy the needs of most extractors defined in the package.
###

### The columns in 'columns' + those involved thru 'filter' and 'orderby' are
### collected and their corresponding tables are INNER JOIN'ed.
TxDb_SELECT_from_INNER_JOIN <- function(txdb, table, columns, filter=list(),
                                        orderby=character(0))
{
    schema_version <- TxDb_schema_version(txdb)
    tables <- TXDB_column2table(columns, from_table=table,
                                schema_version=schema_version)
    where_columns <- names(filter)
    where_tables <- TXDB_column2table(where_columns, from_table=table,
                                      schema_version=schema_version)
    joins <- .TXDB_join_tables(c(table, tables, where_tables))
    orderby_tables <- TXDB_column2table(orderby, from_table=table,
                                        schema_version=schema_version)
    stopifnot(all(orderby_tables %in% .tables_in_joins(joins)))
    use_joins <- length(joins) > 1L
    if (use_joins) {
        columns <- .as_qualified(tables, columns)
        names(filter) <- .as_qualified(where_tables, where_columns)
        orderby <- .as_qualified(orderby_tables, orderby)
    }
    ## .build_SQL_SELECT() uses INNER joins.
    SQL <- .build_SQL_SELECT(columns, joins, distinct=use_joins,
                             filter=filter, orderby=orderby)
    queryAnnotationDb(txdb, SQL)
}

### Can only involve columns (thru 'columns', 'filter', and 'orderby') that
### belong to the tables in TXDB_SPLICING_BUNDLE at the moment.
TxDb_SELECT_from_splicing_bundle <- function(txdb, columns,
                                             filter=list(),
                                             orderby=character(0),
                                             cds_join_type="LEFT")
{
    schema_version <- TxDb_schema_version(txdb)
    tables <- TXDB_column2table(columns, from_table="splicing",
                                schema_version=schema_version)
    where_columns <- names(filter)
    where_tables <- TXDB_column2table(where_columns, from_table="splicing",
                                      schema_version=schema_version)
    orderby_tables <- TXDB_column2table(orderby, from_table="splicing",
                                        schema_version=schema_version)
    joins <- .TXDB_join_splicing_Rtables(c(tables, where_tables,
                                           orderby_tables))
    use_joins <- length(joins) > 1L
    if (use_joins) {
        columns <- .as_qualified(tables, columns)
        names(filter) <- .as_qualified(where_tables, where_columns)
        orderby <- .as_qualified(orderby_tables, orderby)
    }
    from <- .build_SQL_FROM_splicing(joins, cds_join_type=cds_join_type)
    SQL <- .build_SQL_SELECT(columns, from, distinct=FALSE,
                             filter=filter, orderby=orderby)
    queryAnnotationDb(txdb, SQL)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Convenience wrappers to the above flexible helpers for SELECT'ing stuff
### from a given TxDb table
###

TxDb_SELECT_from_chrominfo <- function(txdb, filter=list(),
                                       orderby="_chrom_id")
{
    schema_version <- TxDb_schema_version(txdb)
    columns <- TXDB_table_columns("chrominfo", schema_version=schema_version)
    TxDb_SELECT_from_INNER_JOIN(txdb, "chrominfo", columns,
                                filter=filter, orderby=orderby)
}

TxDb_SELECT_from_transcript <- function(txdb, filter=list(),
                                        orderby="_tx_id")
{
    schema_version <- TxDb_schema_version(txdb)
    columns <- TXDB_table_columns("transcript", schema_version=schema_version)
    TxDb_SELECT_from_INNER_JOIN(txdb, "transcript", columns,
                                filter=filter, orderby=orderby)
}

### Select rows from the virtual table obtained by joining the "splicing",
### "exon", and "cds" tables together.
TxDb_SELECT_from_splicings <- function(txdb, filter=list(),
                                       orderby=c("_tx_id", "exon_rank"),
                                       cds_join_type="LEFT")
{
    schema_version <- TxDb_schema_version(txdb)
    exon_columns <- TXDB_table_columns("exon", schema_version=schema_version)
    cds_columns <- TXDB_table_columns("cds", schema_version=schema_version)
    cds_columns <- cds_columns[c("id", "name", "start", "end")]
    columns <- unique(c("_tx_id", "exon_rank", exon_columns, cds_columns))
    TxDb_SELECT_from_splicing_bundle(txdb, columns,
                                     filter=filter, orderby=orderby,
                                     cds_join_type=cds_join_type)
}

TxDb_SELECT_from_gene <- function(txdb, filter=list(),
                                  orderby=c("_tx_id", "gene_id"))
{
    schema_version <- TxDb_schema_version(txdb)
    columns <- TXDB_table_columns("gene", schema_version=schema_version)
    TxDb_SELECT_from_INNER_JOIN(txdb, "gene", columns,
                                filter=filter, orderby=orderby)
}

