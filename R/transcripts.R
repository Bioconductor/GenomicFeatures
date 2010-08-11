## convert a named list into an SQL where condition
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


## transcripts function and helpers

.newListBySplit <- function(class, x, f)
{
  IRanges:::newCompressedList(class, unlistData = x, splitFactor = f)
}

.geneCharacterList <- function(txdb, tx_ids)
{
  sqlIDs <- paste("(", paste(tx_ids, collapse=","), ")", sep="")
  sql <- paste("SELECT gene_id, _tx_id FROM gene WHERE _tx_id IN ", sqlIDs,
               sep="")
  ans <- dbEasyQuery(txdbConn(txdb), sql)
  ans[["_tx_id"]] <-
    factor(as.character(ans[["_tx_id"]]), levels=as.character(tx_ids))
  ans <- ans[order(ans[["_tx_id"]]), ,drop=FALSE]
  unname(.newListBySplit("CompressedCharacterList", ans[["gene_id"]], ans[["_tx_id"]]))
}

.exonORcdsIntegerList <- function(txdb, tx_ids, type=c("exon", "cds"))
{
  type <- match.arg(type)
  type_id <- paste("_", type, "_id", sep="")
  sqlIDs <- paste("(", paste(tx_ids, collapse=","), ")", sep="")
  sql <- paste("SELECT _tx_id, exon_rank, ", type_id, " ",
               "FROM splicing WHERE _tx_id IN ", sqlIDs, sep="")
  ans <- dbEasyQuery(txdbConn(txdb), sql)
  ans[["_tx_id"]] <-
    factor(as.character(ans[["_tx_id"]]), levels=as.character(tx_ids))
  ans <- ans[order(ans[["_tx_id"]], ans[["exon_rank"]]), ,drop=FALSE]
  .newListBySplit("CompressedIntegerList", ans[[type_id]], ans[["_tx_id"]])
}


transcripts <- function(txdb, vals=NULL, columns=c("tx_id", "tx_name"))
{
  ## check to see if user wanted deprecated function
  if(is.data.frame(txdb))
    stop("Please use 'transcripts_deprecated' for older data.frame-based transcript metadata.")

  ## check that txdb is a TranscriptDb object
  if(!is(txdb,"TranscriptDb"))
    stop("'txdb' must be a TranscriptDb object")

  ## check the vals argument
  validValNames <- c("gene_id", "tx_id", "tx_name", "tx_chrom", "tx_strand")
  if(!is.null(vals) &&
     (!is.list(vals) || is.null(names(vals)) ||
      !all(names(vals) %in% validValNames))) {
    stop("'vals' must be NULL or a list with names being a combination of ",
         paste(dQuote(validValNames), collapse = ", "))
  }
  whichId <- which(names(vals) == "tx_id")
  if(length(whichId) > 0) {
    names(vals)[whichId] <- "transcript._tx_id"
  }

  ## check the columns argument
  validColumns <- c("tx_id", "tx_name", "gene_id", "exon_id","cds_id")
  if(length(columns) > 0 &&
     (!is.character(columns) || !all(columns %in% validColumns))) {
    stop("'columns' must be NULL or a combination of ",
         paste(dQuote(validColumns), collapse = ", "))
  }

  ## create SQL query
  if ("tx_name" %in% columns)
    optionalColumn <- ", tx_name"
  else
    optionalColumn <- ""
  if ("gene_id" %in% names(vals))
    optionalLeftJoin <- "LEFT JOIN gene ON transcript._tx_id=gene._tx_id"
  else
    optionalLeftJoin <- ""
  sql <- paste("SELECT tx_chrom, tx_start, tx_end, tx_strand,",
               "transcript._tx_id AS tx_id", optionalColumn,
               "FROM transcript", optionalLeftJoin,
               .sqlWhereIn(vals),
               "ORDER BY tx_chrom, tx_strand, tx_start, tx_end")

  ## get the data from the database
  ans <- dbEasyQuery(txdbConn(txdb), sql)
  seqlengths <- seqlengths(txdb)
  ans <-
    GRanges(seqnames = factor(ans[["tx_chrom"]], levels = names(seqlengths)),
            ranges = IRanges(start = ans[["tx_start"]],
                             end = ans[["tx_end"]]),
            strand = strand(ans[["tx_strand"]]),
            ans[-c(1:4)],
            seqlengths = seqlengths)

  if(length(ans) > 0 && any(c("gene_id", "exon_id","cds_id") %in% columns)) {
    if("gene_id" %in% columns) {
      elementMetadata(ans)[["gene_id"]] <-
        .geneCharacterList(txdb, elementMetadata(ans)[["tx_id"]])
    }
    if("exon_id" %in% columns) {
      elementMetadata(ans)[["exon_id"]] <-
        .exonORcdsIntegerList(txdb, elementMetadata(ans)[["tx_id"]], "exon")
    }
    if("cds_id" %in% columns) {
      elementMetadata(ans)[["cds_id"]] <-
        .exonORcdsIntegerList(txdb, elementMetadata(ans)[["tx_id"]], "cds")
    }
  }

  if (!("tx_id" %in% columns))
    elementMetadata(ans)[["tx_id"]] <- NULL

  ans
}


## exon and cds functions and helper

.exonORcdsGRanges <- function(txdb, vals=NULL, type=c("exon", "cds"))
{
  type <- match.arg(type)

  if (type == "exon" && is.data.frame(txdb))
    stop("Please use 'exons_deprecated' for older data.frame-based transcript metadata.")

  ## check that txdb is a TranscriptDb object
  if(!is(txdb,"TranscriptDb"))
    stop("'txdb' must be a TranscriptDb object")

  ## check the vals argument
  validValNames <- gsub("TYPE", type, c("TYPE_id", "TYPE_chrom", "TYPE_strand"))
  if(!is.null(vals) &&
     (!is.list(vals) || is.null(names(vals)) ||
      !all(names(vals) %in% validValNames))) {
    stop("'vals' must be NULL or a list with names being a combination of ",
         paste(dQuote(validValNames), collapse = ", "))
  }
  whichId <- which(names(vals) == paste(type, "_id", sep=""))
  if(length(whichId) > 0) {
      names(vals)[whichId] <- paste("TYPE._", type, "_id", sep="")
  }

  ## create base SQL query
  sql <- paste("SELECT TYPE_chrom, TYPE_start, TYPE_end, TYPE_strand,",
               "TYPE._TYPE_id AS TYPE_id FROM TYPE",
               .sqlWhereIn(vals),
               "ORDER BY TYPE_chrom, TYPE_strand, TYPE_start, TYPE_end")
  sql <- gsub("TYPE", type, sql)

  ## get the data from the database
  ans <- dbEasyQuery(txdbConn(txdb), sql)
  seqlengths <- seqlengths(txdb)
  ans <-
    GRanges(seqnames =
            factor(ans[[paste(type, "_chrom", sep="")]],
                   levels = names(seqlengths)),
            ranges = IRanges(start = ans[[paste(type, "_start", sep="")]],
                             end = ans[[paste(type, "_end", sep="")]]),
            strand = strand(ans[[paste(type, "_strand", sep="")]]),
            "TYPE_id" = ans[[paste(type, "_id", sep="")]],
            seqlengths = seqlengths)
  colnames(elementMetadata(ans)) <-
    gsub("TYPE", type, colnames(elementMetadata(ans)))

  ans
}

exons.OLD <- function(txdb, vals=NULL)
{
  .exonORcdsGRanges(txdb, vals=vals, type="exon")
}


cds <- function(txdb, vals=NULL)
{
  .exonORcdsGRanges(txdb, vals=vals, type="cds")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "exons" function.
###
### If we look at the db from an exon centric point of view, then the graph
### of relations between the tables becomes a tree where the 'exon' table is
### the root:
###                               exon
###                                |
###                             splicing
###                            /   |    \
###                   transcript  gene  cds
###

### Each col defined in the db is assigned to the closest table (starting
### from the 'exon' table) where it can be found. The 'exon' element must be
### the 1st in the list:
.EXON.COLMAP <- list(
    exon=c("exon_id", makeFeatureColnames("exon")),
    splicing=c("tx_id", "_tx_id", "exon_rank", "cds_id", "_cds_id"),
    transcript=c("tx_name", "tx_chrom", "tx_strand", "tx_start", "tx_end"),
    gene="gene_id",
    cds=c("cds_name", "cds_chrom", "cds_strand", "cds_start", "cds_end")
)

### For each table that is not the root or a leaf table in the above tree,
### we list the tables that are below it:
.EXON.CHILDTABLES <- list(
    splicing=c("transcript", "gene", "cds")
)

### For each table that is not the root table in the above tree, we specify
### the join condition with the parent table:
.EXON.JOINS <- c(
    splicing="exon._exon_id=splicing._exon_id",
    transcript="splicing._tx_id=transcript._tx_id",
    gene="splicing._tx_id=gene._tx_id",
    cds="splicing._cds_id=cds._cds_id"
)

.exon.getClosestTable <- function(colnames)
{
    ans <- character(length(colnames))
    ans[] <- NA_character_
    for (tablename in names(.EXON.COLMAP))
        ans[colnames %in% .EXON.COLMAP[[tablename]]] <- tablename
    ans
}

.exon.assignColToClosestTable <- function(colnames)
{
    lapply(.EXON.COLMAP, function(cols) intersect(cols, colnames))
}

.exon.asQualifiedColnames <- function(colnames)
{
    colnames[colnames == "tx_id"] <- "_tx_id"
    colnames[colnames == "exon_id"] <- "_exon_id"
    colnames[colnames == "cds_id"] <- "_cds_id"
    paste(.exon.getClosestTable(colnames), colnames, sep=".")
}

.exon.makeSQLfrom <- function(right_tables)
{
    all_tables <- names(.EXON.COLMAP)
    ans <- all_tables[1L]
    for (i in seq_len(length(all_tables))[-1L]) {
        tablename <- all_tables[i]
        children <- c(tablename, .EXON.CHILDTABLES[[tablename]])
        if (length(intersect(right_tables, children)) != 0L)
            ans <- paste(ans, "LEFT JOIN", tablename,
                              "ON", .EXON.JOINS[[tablename]])
    }
    ans
}

.exon.makeSQL <- function(what_cols, right_tables, vals, orderby_cols=NULL)
{
    SQL_what <- paste(.exon.asQualifiedColnames(what_cols), collapse=", ")
    SQL_from <- .exon.makeSQLfrom(right_tables)
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

.exon.checkargColumns <- function(columns)
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

.exon.getWhereCols <- function(vals)
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

.exon.extractPrimaryData <- function(txdb, vals, primary_cols)
{
    CORE.COLUMNS <- c("exon_id", "exon_chrom", "exon_strand",
                      "exon_start", "exon_end")
    what_cols <- unique(c(CORE.COLUMNS, primary_cols))
    where_cols <- .exon.getWhereCols(vals)
    if (is.list(vals))
        names(vals) <- .exon.asQualifiedColnames(where_cols)
    where_tables <- unique(.exon.getClosestTable(where_cols))
    orderby_cols <- c("exon_chrom", "exon_strand", "exon_start", "exon_end")
    .exon.makeSQL(what_cols, where_tables, vals, orderby_cols)
}

.exon.extractForeignData <- function(txdb, ids, assigned_cols)
{
    ans <- NULL
    all_tables <- names(assigned_cols)
    for (i in seq_len(length(all_tables))[-1L]) {
        foreign_cols <- assigned_cols[[i]]
        if (length(foreign_cols) == 0L)
            next
        right_table <- all_tables[i]
        what_cols <- c("exon_id", foreign_cols)
        vals <- list(exon_id=ids)
        names(vals) <- .exon.asQualifiedColnames(names(vals))
        data0 <- .exon.makeSQL(what_cols, right_table, vals)
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

### Note that the current naming of the args is a little bit confusing
### because we have the 'vals' and 'columns' args and it is the latter that is
### related to the "values" slot of the returned GRanges object, not the
### former. More precisely, the user specifies what goes into the "values"
### slot (aka "elementMetadata" slot) thru the 'columns' arg and not thru the
### 'vals' arg.
### TODO:
###   - Rename the 'vals' arg -> 'filter'.
###   - Rename the 'columns' arg -> 'colnames'.
exons <- function(txdb, vals=NULL, columns="exon_id")
#exons <- function(txdb, vals=NULL, columns=c("exon_id", "exon_name"))
{
    if (is.data.frame(txdb))
        stop("Please use 'exons_deprecated' for older data.frame-based ",
             "transcript metadata.")
    if (!is(txdb, "TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")
    .exon.checkargColumns(columns)
    assigned_cols <- .exon.assignColToClosestTable(columns)
    ## Extract the data from the db.
    primary_data <- .exon.extractPrimaryData(txdb, vals, assigned_cols$exon)
    foreign_data <- .exon.extractForeignData(txdb,
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
    ans_values <- c(DataFrame(primary_data[assigned_cols$exon]), foreign_data)
    values(ans) <- ans_values[columns]
    ans
}

