### =========================================================================
### TxDb schema
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


TXDB_CHROMINFO_COLUMNS <- c(
    "_chrom_id",
    "chrom",
    "length",
    "is_circular"
)

TXDB_TRANSCRIPT_COLUMNS <- c(
    id="_tx_id",
    name="tx_name",
    type="tx_type",
    chrom="tx_chrom",
    strand="tx_strand",
    start="tx_start",
    end="tx_end"
)

TXDB_EXON_COLUMNS <- c(
    id="_exon_id",
    name="exon_name",
    chrom="exon_chrom",
    strand="exon_strand",
    start="exon_start",
    end="exon_end"
)

TXDB_CDS_COLUMNS <- c(
    id="_cds_id",
    name="cds_name",
    chrom="cds_chrom",
    strand="cds_strand",
    start="cds_start",
    end="cds_end"
)

TXDB_SPLICING_COLUMNS <- c(
    "_tx_id",
    "exon_rank",
    "_exon_id",
    "_cds_id"
)

TXDB_GENE_COLUMNS <- c(
    "gene_id",
    "_tx_id"
)

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

### Tables "transcript", "exon", and "cds", must have these tags (at a minimum).
TXDB_CORE_TAGS <- c("id", "chrom", "strand", "start", "end")

### The "splicing right tables" can be bundled to the "splicing" table with
### a LEFT JOIN using the TXDB_SPLICING_JOIN_USING columns.
TXDB_SPLICING_RTABLES <- c("transcript", "exon", "cds")
TXDB_SPLICING_JOIN_USING <- setNames(c("_tx_id", "_exon_id", "_cds_id"),
                                     TXDB_SPLICING_RTABLES)
TXDB_SPLICING_BUNDLE <- c("splicing", TXDB_SPLICING_RTABLES)

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

