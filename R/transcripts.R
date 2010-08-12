### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A high-level representation of the relational organization of the db.
###

### THE TRANSCRIPT CENTRIC POINT OF VIEW:
### If we look at the db from a transcript centric point of view, then the
### graph of relations between the tables becomes a tree where the 'transcript'
### table is the root:
###                            transcript
###                             /      \
###                          gene    splicing
###                                  /      \
###                                exon     cds
###
.TRANSCRIPT_CENTRIC_DBDESC <- list(
    CORECOLS=paste("tx", c("id", "chrom", "strand", "start", "end"),
                   sep="_"),
    ### Each col defined in the db is assigned to the closest table (starting
    ### from the 'exon' table) where it can be found. The 'exon' element must
    ### be the 1st in the list:
    COLMAP=list(
        transcript=c("tx_id", makeFeatureColnames("tx")),
        gene="gene_id",
        splicing=c("exon_rank", "exon_id", "_exon_id", "cds_id", "_cds_id"),
        exon=c("exon_name", "exon_chrom", "exon_strand",
               "exon_start", "exon_end"),
        cds=c("cds_name", "cds_chrom", "cds_strand", "cds_start", "cds_end")
    ),
    ### For each table that is not the root or a leaf table in the above tree,
    ### we list the tables that are below it:
    CHILDTABLES=list(
        splicing=c("exon", "cds")
    ),
    ### For each table that is not the root table in the above tree, we specify
    ### the join condition with the parent table:
    JOINS=c(
        gene="transcript._tx_id=gene._tx_id",
        splicing="transcript._tx_id=splicing._tx_id",
        exon="splicing._exon_id=exon._exon_id",
        cds="splicing._cds_id=cds._cds_id"
    )
)

### THE EXON CENTRIC POINT OF VIEW:
### If we look at the db from an exon centric point of view, then the
### graph of relations between the tables becomes a tree where the 'exon'
### table is the root:
###                               exon
###                                |
###                             splicing
###                            /   |    \
###                   transcript  gene  cds
###
.EXON_CENTRIC_DBDESC <- list(
    CORECOLS=paste("exon", c("id", "chrom", "strand", "start", "end"),
                   sep="_"),
    ### Each col defined in the db is assigned to the closest table (starting
    ### from the 'exon' table) where it can be found. The 'exon' element must
    ### be the 1st in the list:
    COLMAP=list(
        exon=c("exon_id", makeFeatureColnames("exon")),
        splicing=c("tx_id", "_tx_id", "exon_rank", "cds_id", "_cds_id"),
        transcript=c("tx_name", "tx_chrom", "tx_strand", "tx_start", "tx_end"),
        gene="gene_id",
        cds=c("cds_name", "cds_chrom", "cds_strand", "cds_start", "cds_end")
    ),
    ### For each table that is not the root or a leaf table in the above tree,
    ### we list the tables that are below it:
    CHILDTABLES=list(
        splicing=c("transcript", "gene", "cds")
    ),
    ### For each table that is not the root table in the above tree, we specify
    ### the join condition with the parent table:
    JOINS=c(
        splicing="exon._exon_id=splicing._exon_id",
        transcript="splicing._tx_id=transcript._tx_id",
        gene="splicing._tx_id=gene._tx_id",
        cds="splicing._cds_id=cds._cds_id"
    )
)

### THE CDS CENTRIC POINT OF VIEW:
### If we look at the db from a cds centric point of view, then the
### graph of relations between the tables becomes a tree where the 'cds'
### table is the root:
###                               cds
###                                |
###                             splicing
###                            /   |    \
###                   transcript  gene  exon
###
.CDS_CENTRIC_DBDESC <- list(
    CORECOLS=paste("cds", c("id", "chrom", "strand", "start", "end"),
                   sep="_"),
    ### Each col defined in the db is assigned to the closest table (starting
    ### from the 'exon' table) where it can be found. The 'exon' element must
    ### be the 1st in the list:
    COLMAP=list(
        cds=c("cds_id", makeFeatureColnames("cds")),
        splicing=c("tx_id", "_tx_id", "exon_rank", "exon_id", "_exon_id"),
        transcript=c("tx_name", "tx_chrom", "tx_strand", "tx_start", "tx_end"),
        gene="gene_id",
        exon=c("exon_name", "exon_chrom", "exon_strand",
               "exon_start", "exon_end")
    ),
    ### For each table that is not the root or a leaf table in the above tree,
    ### we list the tables that are below it:
    CHILDTABLES=list(
        splicing=c("transcript", "gene", "exon")
    ),
    ### For each table that is not the root table in the above tree, we specify
    ### the join condition with the parent table:
    JOINS=c(
        splicing="cds._cds_id=splicing._cds_id",
        transcript="splicing._tx_id=transcript._tx_id",
        gene="splicing._tx_id=gene._tx_id",
        exon="splicing._exon_id=exon._exon_id"
    )
)

.DBDESC <- list(
    transcript=.TRANSCRIPT_CENTRIC_DBDESC,
    exon=.EXON_CENTRIC_DBDESC,
    cds=.CDS_CENTRIC_DBDESC
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utility functions for accessing the high-level representation of the
### relational organization of the db.
###

.getClosestTable <- function(root_table, colnames)
{
    COLMAP <- .DBDESC[[root_table]]$COLMAP
    ans <- character(length(colnames))
    ans[] <- NA_character_
    for (tablename in names(COLMAP))
        ans[colnames %in% COLMAP[[tablename]]] <- tablename
    ans
}

.assignColToClosestTable <- function(root_table, colnames)
{
    COLMAP <- .DBDESC[[root_table]]$COLMAP
    lapply(COLMAP, function(cols) intersect(cols, colnames))
}

.asQualifiedColnames <- function(root_table, colnames)
{
    colnames[colnames == "tx_id"] <- "_tx_id"
    colnames[colnames == "exon_id"] <- "_exon_id"
    colnames[colnames == "cds_id"] <- "_cds_id"
    paste(.getClosestTable(root_table, colnames), colnames, sep=".")
}

.makeSQLfrom <- function(root_table, right_tables)
{
    COLMAP <- .DBDESC[[root_table]]$COLMAP
    CHILDTABLES <- .DBDESC[[root_table]]$CHILDTABLES
    JOINS <- .DBDESC[[root_table]]$JOINS
    all_tables <- names(COLMAP)
    ans <- all_tables[1L]
    for (i in seq_len(length(all_tables))[-1L]) {
        tablename <- all_tables[i]
        children <- c(tablename, CHILDTABLES[[tablename]])
        if (length(intersect(right_tables, children)) != 0L)
            ans <- paste(ans, "LEFT JOIN", tablename,
                              "ON", JOINS[[tablename]])
    }
    ans
}

### convert a named list into an SQL where condition
.sqlWhereIn <- function(vals)
{
    if (length(vals) == 0L)
        return("")
    sql <-
      lapply(seq_len(length(vals)), function(i) {
               v <- vals[[i]]
               if (!is.numeric(v))
                 v <- paste("'", v, "'", sep="")
               v <- paste("(", paste(v, collapse=","), ")", sep="")
               v <- paste(names(vals)[i], " IN ", v, sep="")
               paste("(", v, ")", sep="")
            })
    paste("WHERE", paste(unlist(sql), collapse = " AND "))
}

.extractData <- function(root_table, txdb, what_cols, right_tables, vals,
                         orderby_cols=NULL)
{
    SQL_what <- paste(.asQualifiedColnames(root_table, what_cols),
                      collapse=", ")
    SQL_from <- .makeSQLfrom(root_table, right_tables)
    SQL_where <- .sqlWhereIn(vals)
    if (length(orderby_cols) == 0L)
        SQL_orderby <- ""
    else
        SQL_orderby <- paste("ORDER BY", paste(orderby_cols, collapse=", "))
    SQL <- paste("SELECT DISTINCT", SQL_what, "FROM", SQL_from,
                 SQL_where, SQL_orderby)
    ans <- dbEasyQuery(txdbConn(txdb), SQL)
    names(ans) <- what_cols
    ans
}

.getWhereCols <- function(vals)
{
    if (is.null(vals))
        return(character(0))
    ans <- NULL
    if (is.list(vals)) {
        if (length(vals) == 0L)
            return(character(0))
        ans <- names(vals)
    }
    if (is.null(ans))
        stop("'vals' must be NULL or a named list")
    VALID.COLUMNS <- c("gene_id",
                       "tx_id", "tx_name", "tx_chrom", "tx_strand",
                       "exon_id", "exon_name", "exon_chrom", "exon_strand",
                       "cds_id", "cds_name", "cds_chrom", "cds_strand")
    if (!all(ans %in% VALID.COLUMNS)) {
        validColumns <- paste("\"", VALID.COLUMNS, "\"",
                              sep="", collapse = ", ")
        stop("'vals' must be NULL or a list with names ",
             "in ", validColumns)
    }
    ans
}

.extractPrimaryData <- function(root_table, txdb, vals, primary_cols)
{
    CORECOLS <- .DBDESC[[root_table]]$CORECOLS
    orderby_cols <- CORECOLS[2:5]
    what_cols <- unique(c(CORECOLS, primary_cols))
    where_cols <- .getWhereCols(vals)
    if (is.list(vals))
        names(vals) <- .asQualifiedColnames(root_table, where_cols)
    where_tables <- unique(.getClosestTable(root_table, where_cols))
    .extractData(root_table, txdb, what_cols, where_tables, vals, orderby_cols)
}

.extractForeignData <- function(root_table, txdb, ids, assigned_cols)
{
    primary_key <- .DBDESC[[root_table]]$CORECOLS[1L]
    ans <- NULL
    all_tables <- names(assigned_cols)
    for (i in seq_len(length(all_tables))[-1L]) {
        foreign_cols <- assigned_cols[[i]]
        if (length(foreign_cols) == 0L)
            next
        right_table <- all_tables[i]
        what_cols <- c(primary_key, foreign_cols)
        vals <- list(ids)
        names(vals) <- .asQualifiedColnames(root_table, primary_key)
        data0 <- .extractData(root_table, txdb, what_cols, right_table, vals)
        data <- lapply(data0[ , -1L, drop=FALSE],
                       function(col0)
                       {
                           col <- split(col0, data0[[1L]])
                           col <- col[as.character(ids)]
                           class0 <- class(col0)
                           class <- paste(toupper(substr(class0, 1L, 1L)),
                                          substr(class0, 2L, nchar(class0)),
                                          "List", sep="")
                           get(class)(unname(col))
                       })
        data <- DataFrame(data)
        if (is.null(ans))
            ans <- data
        else
            ans <- c(ans, data)
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "transcripts" function.
###

.transcripts.checkargColumns <- function(columns)
{
    if (is.null(columns))
        return()
    VALID.COLUMNS <- c("gene_id",
                       "tx_id", "tx_name",
                       "exon_id", "exon_name", "exon_chrom", "exon_strand",
                       "cds_id", "cds_name", "cds_chrom", "cds_strand")
    if (is.character(columns) && all(columns %in% VALID.COLUMNS))
        return()
    validColumns <- paste("\"", VALID.COLUMNS, "\"", sep="", collapse = ", ")
    stop("'columns' must be NULL or a character vector ",
         "with values in ", validColumns)
}

### Note that the current naming of the args is a little bit confusing
### because we have the 'vals' and 'columns' args and it is the latter that is
### related to the "values" slot of the returned GRanges object, not the
### former. More precisely, the user specifies what goes into the "values"
### slot (aka "elementMetadata" slot) thru the 'columns' arg and not thru the
### 'vals' arg.
### TODO:
###   - Rename the 'vals' arg -> 'filter'.
###   - Rename the 'columns' arg -> 'colnames'.
transcripts <- function(txdb, vals=NULL, columns=c("tx_id", "tx_name"))
{
    if (is.data.frame(txdb))
        stop("Please use 'transcripts_deprecated' for older ",
             "data.frame-based transcript metadata.")
    if (!is(txdb, "TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")
    .transcripts.checkargColumns(columns)
    assigned_cols <- .assignColToClosestTable("transcript", columns)
    ## Extract the data from the db.
    primary_data <- .extractPrimaryData("transcript", txdb,
                        vals, assigned_cols$transcript)
    foreign_data <- .extractForeignData("transcript", txdb,
                        primary_data[["tx_id"]], assigned_cols)
    ## Construct the GRanges object and return it.
    ans_seqlengths <- seqlengths(txdb)
    ans_seqnames <- factor(primary_data[["tx_chrom"]],
                           levels=names(ans_seqlengths))
    ans_ranges <- IRanges(start=primary_data[["tx_start"]],
                          end=primary_data[["tx_end"]])
    ans_strand <- strand(primary_data[["tx_strand"]])
    ans <- GRanges(seqnames=ans_seqnames,
                   ranges=ans_ranges,
                   strand=ans_strand,
                   seqlengths=ans_seqlengths)
    ans_values <- c(DataFrame(primary_data[assigned_cols$transcript]),
                    foreign_data)
    values(ans) <- ans_values[columns]
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "exons" function.
###

.exons.checkargColumns <- function(columns)
{
    if (is.null(columns))
        return()
    VALID.COLUMNS <- c("gene_id",
                       "tx_id", "tx_name", "tx_chrom", "tx_strand",
                       "exon_id", "exon_name",
                       "cds_id", "cds_name", "cds_chrom", "cds_strand")
    if (is.character(columns) && all(columns %in% VALID.COLUMNS))
        return()
    validColumns <- paste("\"", VALID.COLUMNS, "\"", sep="", collapse = ", ")
    stop("'columns' must be NULL or a character vector ",
         "with values in ", validColumns)
}

#exons <- function(txdb, vals=NULL, columns=c("exon_id", "exon_name"))
exons <- function(txdb, vals=NULL, columns="exon_id")
{
    if (is.data.frame(txdb))
        stop("Please use 'exons_deprecated' for older ",
             "data.frame-based transcript metadata.")
    if (!is(txdb, "TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")
    .exons.checkargColumns(columns)
    assigned_cols <- .assignColToClosestTable("exon", columns)
    ## Extract the data from the db.
    primary_data <- .extractPrimaryData("exon", txdb,
                        vals, assigned_cols$exon)
    foreign_data <- .extractForeignData("exon", txdb,
                        primary_data[["exon_id"]], assigned_cols)
    ## Construct the GRanges object and return it.
    ans_seqlengths <- seqlengths(txdb)
    ans_seqnames <- factor(primary_data[["exon_chrom"]],
                           levels=names(ans_seqlengths))
    ans_ranges <- IRanges(start=primary_data[["exon_start"]],
                          end=primary_data[["exon_end"]])
    ans_strand <- strand(primary_data[["exon_strand"]])
    ans <- GRanges(seqnames=ans_seqnames,
                   ranges=ans_ranges,
                   strand=ans_strand,
                   seqlengths=ans_seqlengths)
    ans_values <- c(DataFrame(primary_data[assigned_cols$exon]),
                    foreign_data)
    values(ans) <- ans_values[columns]
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "cds" function.
###

.cds.checkargColumns <- function(columns)
{
    if (is.null(columns))
        return()
    VALID.COLUMNS <- c("gene_id",
                       "tx_id", "tx_name", "tx_chrom", "tx_strand",
                       "exon_id", "exon_name", "exon_chrom", "exon_strand",
                       "cds_id", "cds_name")
    if (is.character(columns) && all(columns %in% VALID.COLUMNS))
        return()
    validColumns <- paste("\"", VALID.COLUMNS, "\"", sep="", collapse = ", ")
    stop("'columns' must be NULL or a character vector ",
         "with values in ", validColumns)
}

#cds <- function(txdb, vals=NULL, columns=c("cds_id", "cds_name"))
cds <- function(txdb, vals=NULL, columns="cds_id")
{
    if (!is(txdb, "TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")
    .cds.checkargColumns(columns)
    assigned_cols <- .assignColToClosestTable("cds", columns)
    ## Extract the data from the db.
    primary_data <- .extractPrimaryData("cds", txdb,
                        vals, assigned_cols$cds)
    foreign_data <- .extractForeignData("cds", txdb,
                        primary_data[["cds_id"]], assigned_cols)
    ## Construct the GRanges object and return it.
    ans_seqlengths <- seqlengths(txdb)
    ans_seqnames <- factor(primary_data[["cds_chrom"]],
                           levels=names(ans_seqlengths))
    ans_ranges <- IRanges(start=primary_data[["cds_start"]],
                          end=primary_data[["cds_end"]])
    ans_strand <- strand(primary_data[["cds_strand"]])
    ans <- GRanges(seqnames=ans_seqnames,
                   ranges=ans_ranges,
                   strand=ans_strand,
                   seqlengths=ans_seqlengths)
    ans_values <- c(DataFrame(primary_data[assigned_cols$cds]),
                    foreign_data)
    values(ans) <- ans_values[columns]
    ans
}

