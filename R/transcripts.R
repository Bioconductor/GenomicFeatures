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

.exon.getClosestTableClosestTable <- function(colnames)
{
    ans <- character(length(colnames))
    ans[] <- NA_character_
    for (tablename in names(.EXON.COLMAP))
        ans[colnames %in% .EXON.COLMAP[[tablename]]] <- tablename
    ans
}

.exon.asQualifiedColnames <- function(colnames)
{
    colnames[colnames == "tx_id"] <- "_tx_id"
    colnames[colnames == "exon_id"] <- "_exon_id"
    colnames[colnames == "cds_id"] <- "_cds_id"
    paste(.exon.getClosestTableClosestTable(colnames), colnames, sep=".")
}

.exon.assignColToClosestTable <- function(colnames)
{
    lapply(.EXON.COLMAP, function(cols) intersect(cols, colnames))
}

.exon.assignExtraColToClosestTable <- function(columns)
{
    if (is.null(columns))
        return(.exon.assignColToClosestTable(character(0)))
    VALID.COLUMNS <- c("gene_id",
                       "tx_id", "tx_name", "tx_chrom", "tx_strand",
                       "exon_id", "exon_name",
                       "cds_id", "cds_name", "cds_chrom", "cds_strand")
    if (!is.character(columns) || !all(columns %in% VALID.COLUMNS)) {
        validColumns <- paste("\"", VALID.COLUMNS, "\"",
                              sep="", collapse = ", ")
        stop("'columns' must be NULL or a character vector ",
             "with values in ", validColumns)
    }
    .exon.assignColToClosestTable(columns)
}

.exon.assignFilterColToClosestTable <- function(columns)
{
    if (is.null(columns))
        return(.exon.assignColToClosestTable(character(0)))
    VALID.COLUMNS <- c("gene_id",
                       "tx_id", "tx_name", "tx_chrom", "tx_strand",
                       "exon_id", "exon_name", "exon_chrom", "exon_strand",
                       "cds_id", "cds_name", "cds_chrom", "cds_strand")
    if (!is.character(columns) || !all(columns %in% VALID.COLUMNS)) {
        validColumns <- paste("\"", VALID.COLUMNS, "\"",
                              sep="", collapse = ", ")
        stop("'vals' must be NULL or a list with names ",
             "in ", validColumns)
    }
    .exon.assignColToClosestTable(columns)
}

.exon.getFilterTables <- function(filter_cols)
{
    filter_cols <- .exon.assignFilterColToClosestTable(filter_cols)
    names(filter_cols)[sapply(filter_cols, length) != 0L]
}

.exon.makeSQLfrom <- function(target_tables)
{
    all_tables <- names(.EXON.COLMAP)
    SQL_from <- all_tables[1L]
    for (i in seq_len(length(all_tables))[-1L]) {
        tablename <- all_tables[i]
        children <- c(tablename, .EXON.CHILDTABLES[[tablename]])
        if (length(intersect(target_tables, children)) != 0L)
            SQL_from <- paste(SQL_from,
                              "LEFT JOIN", tablename,
                              "ON", .EXON.JOINS[[tablename]])
    }
    SQL_from
}

.exon.extractPrimaryData <- function(txdb, vals, primary_cols)
{
    CORE.COLUMNS <- c("exon_id", "exon_chrom", "exon_strand",
                      "exon_start", "exon_end")
    what_cols <- unique(c(CORE.COLUMNS, primary_cols))
    SQL_what <- paste(.exon.asQualifiedColnames(what_cols), collapse=", ")
    filter_cols <- NULL
    if (is.null(vals)) {
        filter_cols <- character(0)
    } else if (is.list(vals)) {
        if (length(vals) == 0L)
            filter_cols <- character(0)
        else
            filter_cols <- names(vals)
    }
    if (is.null(filter_cols))
        stop("'vals' must be NULL or a named list")
    if (is.list(vals))
        names(vals) <- .exon.asQualifiedColnames(filter_cols)
    filter_tables <- .exon.getFilterTables(filter_cols)
    SQL_from <- .exon.makeSQLfrom(filter_tables)
    SQL_where <- .sqlWhereIn(vals)
    SQL <- paste("SELECT DISTINCT", SQL_what, "FROM", SQL_from, SQL_where,
                 "ORDER BY exon_chrom, exon_strand, exon_start, exon_end")
    data <- dbEasyQuery(txdbConn(txdb), SQL)
    names(data) <- what_cols
    data
}

.exon.extractForeignData <- function(txdb, ids, extra_cols)
{
    ans <- NULL
    all_tables <- names(extra_cols)
    for (i in seq_len(length(all_tables))[-1L]) {
        foreign_cols <- extra_cols[[i]]
        if (length(foreign_cols) == 0L)
            next
        tablename <- all_tables[i]
        what_cols <- c("exon_id", foreign_cols)
        SQL_what <- paste(.exon.asQualifiedColnames(what_cols), collapse=", ")
        SQL_from <- .exon.makeSQLfrom(tablename)
        vals <- list(exon_id=ids)
        names(vals) <- .exon.asQualifiedColnames(names(vals))
        SQL_where <- .sqlWhereIn(vals)
        SQL <- paste("SELECT DISTINCT", SQL_what, "FROM", SQL_from, SQL_where)
        data0 <- dbEasyQuery(txdbConn(txdb), SQL)
        names(data0) <- what_cols
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
    extra_cols <- .exon.assignExtraColToClosestTable(columns)
    ## Extract the data from the db.
    primary_data <- .exon.extractPrimaryData(txdb, vals, extra_cols$exon)
    foreign_data <- .exon.extractForeignData(txdb,
                        primary_data[["exon_id"]], extra_cols)
    ## Construct the GRanges object to return.
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
    ans_values <- c(DataFrame(primary_data[extra_cols$exon]), foreign_data)
    values(ans) <- ans_values[columns]
    ans
}

