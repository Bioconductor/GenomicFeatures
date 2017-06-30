### =========================================================================
### TxDb schema
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###
### 7 tables:
###   - chrominfo
###   - transcript
###   - exon
###   - cds
###   - splicing
###   - gene
###   - metadata (not described here)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Table columns
###

### 'chrominfo' table

TXDB_CHROMINFO_COLDEFS <- c(
    `_chrom_id`="INTEGER PRIMARY KEY",
    chrom="TEXT UNIQUE NOT NULL",
    length="INTEGER NULL",
    is_circular="INTEGER NULL"
)

TXDB_CHROMINFO_COLUMNS <- names(TXDB_CHROMINFO_COLDEFS)

### 'transcript', 'exon', and 'cds' tables (a.k.a. "feature tables")

TXDB_FEATURE_COLDEFS <- c(
    id="INTEGER PRIMARY KEY",
    name="TEXT NULL",
    type="TEXT NULL",
    chrom="TEXT NOT NULL",
    strand="TEXT NOT NULL",
    start="INTEGER NOT NULL",
    end="INTEGER NOT NULL"
)

### Tables "transcript", "exon", and "cds" must at least have columns with
### the core column tags.
TXDB_CORE_COLTAGS <- c("id", "chrom", "strand", "start", "end")
TXDB_ALL_COLTAGS <- names(TXDB_FEATURE_COLDEFS)
TXDB_EXON_OR_CDS_COLTAGS <- TXDB_ALL_COLTAGS[TXDB_ALL_COLTAGS != "type"]

.make_feature_columns <- function(prefix, tags)
{
    fmt <- paste0("%s_", tags)
    id_pos <- match("id", tags)
    stopifnot(identical(id_pos, 1L))
    fmt[[id_pos]] <- paste0("_", fmt[[id_pos]])
    setNames(sprintf(fmt, prefix), tags)
}

TXDB_TRANSCRIPT_COLUMNS <- .make_feature_columns("tx", TXDB_ALL_COLTAGS)
TXDB_EXON_COLUMNS <- .make_feature_columns("exon", TXDB_EXON_OR_CDS_COLTAGS)
TXDB_CDS_COLUMNS <- .make_feature_columns("cds", TXDB_EXON_OR_CDS_COLTAGS)

### 'splicing' table

TXDB_SPLICING_COLDEFS <- c(
    `_tx_id`="INTEGER NOT NULL",
    exon_rank="INTEGER NOT NULL",
    `_exon_id`="INTEGER NOT NULL",
    `_cds_id`="INTEGER NULL"
)

TXDB_SPLICING_COLUMNS <- names(TXDB_SPLICING_COLDEFS)

### 'gene' table

TXDB_GENE_COLDEFS <- c(
    gene_id="TEXT NOT NULL",
    `_tx_id`="INTEGER NOT NULL"
)

TXDB_GENE_COLUMNS <- names(TXDB_GENE_COLDEFS)


### Order of tables matters! "transcript" must be before "splicing" and "gene",
### and "exon" and "cds" must be before "splicing". See TXDB_column2table()
### below why.
TXDB_COLUMNS <- list(
    chrominfo=TXDB_CHROMINFO_COLUMNS,
    transcript=TXDB_TRANSCRIPT_COLUMNS,
    exon=TXDB_EXON_COLUMNS,
    cds=TXDB_CDS_COLUMNS,
    splicing=TXDB_SPLICING_COLUMNS,
    gene=TXDB_GENE_COLUMNS
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Build CREATE TABLE statements
###

.build_SQL_CREATE_TABLE <- function(table, coldefs, constraints=NULL)
{
    SQL <- "CREATE TABLE %s (%s\n)"
    coldefs <- c(paste(names(coldefs), coldefs), constraints)
    coldefs <- paste("\n  ", coldefs, collapse=",")
    sprintf(SQL, table, coldefs)
}

build_SQL_CREATE_chrominfo_table <- function()
{
    .build_SQL_CREATE_TABLE("chrominfo", TXDB_CHROMINFO_COLDEFS)
}

build_SQL_CREATE_feature_table <- function(table)
{
    columns <- TXDB_COLUMNS[[table]]
    coldefs <- setNames(TXDB_FEATURE_COLDEFS[names(columns)], columns)
    foreign_key <- sprintf("FOREIGN KEY (%s) REFERENCES chrominfo (chrom)",
                           columns[["chrom"]])
    .build_SQL_CREATE_TABLE(table, coldefs, foreign_key)
}

build_SQL_CREATE_splicing_table <- function()
{
    unique_key <- "UNIQUE (_tx_id, exon_rank)"
    foreign_keys <- sprintf("FOREIGN KEY (_%s_id) REFERENCES %s",
                            c("tx", "exon", "cds"),
                            c("transcript", "exon", "cds"))
    constraints <- c(unique_key, foreign_keys)
    .build_SQL_CREATE_TABLE("splicing", TXDB_SPLICING_COLDEFS, constraints)
}

build_SQL_CREATE_gene_table <- function()
{
    unique_key <- "UNIQUE (gene_id, _tx_id)"
    foreign_key <- "FOREIGN KEY (_tx_id) REFERENCES transcript"
    constraints <- c(unique_key, foreign_key)
    .build_SQL_CREATE_TABLE("gene", TXDB_GENE_COLDEFS, constraints)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Relationship between the 'splicing' table and the "feature tables"
###
### The 'splicing' table is the glue between the "feature tables".
###

### The "splicing right tables" can be bundled to the "splicing" table with
### a LEFT JOIN using the TXDB_SPLICING_JOIN_USING columns.
TXDB_SPLICING_RTABLES <- c("transcript", "exon", "cds")
TXDB_SPLICING_JOIN_USING <- setNames(c("_tx_id", "_exon_id", "_cds_id"),
                                     TXDB_SPLICING_RTABLES)
TXDB_SPLICING_BUNDLE <- c("splicing", TXDB_SPLICING_RTABLES)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions
###

TXDB_tables <- function() names(TXDB_COLUMNS)

TXDB_table_columns <- function(table, schema_version=NA)
{
    columns <- TXDB_COLUMNS[[table]]
    if (is.na(schema_version))
        return(columns)
    if (table == "transcript" && schema_version < numeric_version("1.1"))
        columns <- columns[columns != "tx_type"]
    columns
}

### When the same column belongs to more than one table (e.g. "_tx_id",
### "_exon_id", or "_cds_id"), then the table for which the column is a
### primary key is chosen by default. This behavior can be changed by passing
### the name of a table to 'from_table' in which case the priority is given to
### that table.
TXDB_column2table <- function(columns, from_table=NA, schema_version=NA)
{
    if (length(columns) == 0L)
        return(character(0))
    tables <- sapply(columns,
        function(column) {
            for (table in TXDB_tables()) {
                table_columns <- TXDB_table_columns(table,
                                     schema_version=schema_version)
                if (column %in% table_columns)
                    return(table)
            }
            if (is.na(schema_version)) {
                in_schema <- ""
            } else {
                in_schema <- c(" in db schema ", as.character(schema_version))
            }
            stop(column, ": no such column", in_schema)
        }
    )
    if (!is.na(from_table)) {
        table_columns <- TXDB_table_columns(from_table,
                                            schema_version=schema_version)
        tables[columns %in% table_columns] <- from_table
    }
    tables
}

