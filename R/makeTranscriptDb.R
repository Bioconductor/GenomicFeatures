### =========================================================================
### Making TranscriptDb objects
### -------------------------------------------------------------------------


.DB_TYPE_NAME <- "Db type"
.DB_TYPE_VALUE <- "TranscriptDb"  # same as the name of the class

.makeFeatureColnames <- function(feature_shortname)
{
    suffixes <- c("_id", "_name", "_chrom", "_strand", "_start", "_end")
    prefixes <- c("_", rep.int("", length(suffixes) - 1L))
    paste(prefixes, feature_shortname, suffixes, sep="")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity of a TranscriptDb object.
###

### The specified table must have *at least* the cols specified in 'colnames'.
### It's OK if it has more cols or if it has them in a different order.
.valid.table.colnames <- function(conn, tablename, colnames)
{
    tmp <- try(dbExistsTable(conn, tablename), silent=TRUE)
    if (is(tmp, "try-error"))
        return("invalid DB file")
    if (!tmp)
        return(paste("the DB has no ", tablename, " table", sep=""))
    sql0 <- paste("SELECT * FROM ", tablename, " LIMIT 0", sep="")
    data0 <- dbGetQuery(conn, sql0)
    colnames0 <- colnames(data0)
    if (!all(colnames %in% colnames0)) {
        msg <- paste("the ", tablename, " table in the DB doesn't have ",
                     "all the expected columns (",
                     paste("\"", colnames, "\"", sep="", collapse=", "),
                     ")", sep="")
        return(msg)
    }
    NULL
}

.valid.metadata.table <- function(conn)
{
    colnames <- c("name", "value")
    msg <- .valid.table.colnames(conn, "metadata", colnames)
    if (!is.null(msg))
        return(msg)
    sql <- paste("SELECT * FROM metadata",
                 " WHERE name = '", .DB_TYPE_NAME, "'", sep="")
    data <- dbGetQuery(conn, sql)
    if (nrow(data) != 1L) {
        msg <- paste("the metadata table in the DB has 0 or more ",
                     "than 1 '", .DB_TYPE_NAME, "' entries", sep="")
        return(msg)
    }
    db_type <- data[["value"]]
    if (is.na(db_type) || db_type != .DB_TYPE_VALUE) {
        msg <- paste("'", .DB_TYPE_NAME, "' is not \"", .DB_TYPE_VALUE,
                     "\"", sep="")
        return(msg)
    }
    NULL
}

### TODO: Add more checks!
.valid.transcript.table <- function(conn)
{
    colnames <- .makeFeatureColnames("tx")
    msg <- .valid.table.colnames(conn, "transcript", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

### TODO: Add more checks!
.valid.exon.table <- function(conn)
{
    colnames <- .makeFeatureColnames("exon")
    msg <- .valid.table.colnames(conn, "exon", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

### TODO: Add more checks!
.valid.cds.table <- function(conn)
{
    colnames <- .makeFeatureColnames("cds")
    msg <- .valid.table.colnames(conn, "cds", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

### TODO: Add more checks!
.valid.splicing.table <- function(conn)
{
    colnames <- c("_tx_id", "exon_rank", "_exon_id", "_cds_id")
    msg <- .valid.table.colnames(conn, "splicing", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

### TODO: Add more checks!
.valid.gene.table <- function(conn)
{
    colnames <- c("gene_id", "_tx_id")
    msg <- .valid.table.colnames(conn, "gene", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

### TODO: Add more checks!
.valid.chrominfo.table <- function(conn)
{
    colnames <- c("_chrom_id", "chrom", "length")
    msg <- .valid.table.colnames(conn, "chrominfo", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

.valid.TranscriptDb <- function(x)
{
    c(.valid.metadata.table(x@conn),
      .valid.transcript.table(x@conn),
      .valid.exon.table(x@conn),
      .valid.cds.table(x@conn),
      .valid.splicing.table(x@conn),
      .valid.gene.table(x@conn),
      .valid.chrominfo.table(x@conn))
}

setValidity2("TranscriptDb", .valid.TranscriptDb)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Saving/loading.
###

saveFeatures <- function(x, file)
{
    if (!is(x, "TranscriptDb"))
        stop("'x' must be a TranscriptDb object")
    if (!isSingleString(file))
        stop("'file' must be a single string")
    sqliteCopyDatabase(x@conn, file)
}

loadFeatures <- function(file)
{
    if (!isSingleString(file))
        stop("'file' must be a single string")
    if(!file.exists(file))
        stop("file '", file, "' does not exist")
    conn <- dbConnect(SQLite(), file)
    new("TranscriptDb", conn=conn)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions for makeTranscriptDb().
###

.checkargColnames <- function(arg, required_colnames, optional_colnames,
                              argname)
{
    supported_colnames <- c(required_colnames, optional_colnames)
    if (!is.data.frame(arg))
        stop("'", argname, "' must be a data frame")
    if (!all(required_colnames %in% colnames(arg)))
        stop("'", argname, "' must have at least the following cols: ",
             paste(required_colnames, collapse=", "))
    is_supported_col <- colnames(arg) %in% supported_colnames
    if (any(duplicated(colnames(arg)[is_supported_col])))
        stop("'", argname, "' has duplicated colnames")
    if (any(!is_supported_col))
        warning("ignoring the following cols in '", argname, "': ",
            paste(colnames(arg)[!is_supported_col], collapse=", "))
}

.isCharacterVectorOrFactor <- function(x)
{
    is.character(x) || (is.factor(x) && is.character(levels(x)))
}

.checkForeignKey <- function(referring_vals, referring_type, referring_colname,
                             referred_vals, referred_type, referred_colname)
{
    if (!is.na(referring_type) && !is(referring_vals, referring_type))
        stop("'", referring_colname, "' must be of type ", referring_type)
    if (!is.na(referred_type) && !is(referred_vals, referred_type))
        stop("'", referred_colname, "' must be of type ", referred_type)
    if (any(is.na(referring_vals)))
        stop("'", referring_colname, "' cannot contain NAs")
    if (!all(referring_vals %in% referred_vals))
        stop("all the values in '", referring_colname, "' must ",
             "be present in '", referred_colname, "'")
}

.normargTranscripts <- function(transcripts)
{
    .REQUIRED_COLS <- c("tx_id", "tx_chrom", "tx_strand", "tx_start", "tx_end")
    .OPTIONAL_COLS <- "tx_name"
    .checkargColnames(transcripts, .REQUIRED_COLS, .OPTIONAL_COLS,
                      "transcripts")
    ## Check 'tx_id'.
    if (!is.integer(transcripts$tx_id) || any(is.na(transcripts$tx_id)))
        stop("'transcripts$tx_id' must be an integer vector, with no NAs")
    if (any(duplicated(transcripts$tx_id)))
        stop("'transcripts$tx_id' contains duplicated values")
    ## Check 'tx_name'.
    if ("tx_name" %in% colnames(transcripts)
     && !.isCharacterVectorOrFactor(transcripts$tx_name))
        stop("'transcripts$tx_name' must be a character vector (or factor)")
    ## Check 'tx_chrom'.
    if (!.isCharacterVectorOrFactor(transcripts$tx_chrom)
     || any(is.na(transcripts$tx_chrom)))
        stop("'transcripts$tx_chrom' must be a character vector (or factor) ",
             "with no NAs")
    ## Check 'tx_strand'.
    if (!.isCharacterVectorOrFactor(transcripts$tx_strand)
     || any(is.na(transcripts$tx_strand)))
        stop("'transcripts$tx_strand' must be a character vector (or factor) ",
             "with no NAs")
    if (!all(transcripts$tx_strand %in% c("+", "-")))
        stop("values in 'transcripts$tx_strand' must be \"+\" or \"-\"")
    ## Check 'tx_start'.
    if (!is.numeric(transcripts$tx_start)
     || any(is.na(transcripts$tx_start)))
        stop("'transcripts$tx_start' must be an integer vector with no NAs")
    if (!is.integer(transcripts$tx_start))
        transcripts$tx_start <- as.integer(transcripts$tx_start)
    ## Check 'tx_end'.
    if (!is.numeric(transcripts$tx_end)
     || any(is.na(transcripts$tx_end)))
        stop("'transcripts$tx_end' must be an integer vector with no NAs")
    if (!is.integer(transcripts$tx_end))
        transcripts$tx_end <- as.integer(transcripts$tx_end)
    ## Check 'tx_start <= tx_end'.
    if (any(transcripts$tx_start > transcripts$tx_end))
        stop("transcript starts must be <= transcript ends")
    transcripts
}

.normargSplicings <- function(splicings, unique_tx_ids)
{
    .REQUIRED_COLS <- c("tx_id", "exon_rank", "exon_start", "exon_end")
    .OPTIONAL_COLS <- c("exon_id", "exon_name", "exon_chrom", "exon_strand",
                        "cds_id", "cds_name", "cds_start", "cds_end")
    .checkargColnames(splicings, .REQUIRED_COLS, .OPTIONAL_COLS, "splicings")
    ## Check 'tx_id'.
    .checkForeignKey(splicings$tx_id, "integer", "splicings$tx_id",
                     unique_tx_ids, "integer", "transcripts$tx_id")
    ## Check 'exon_rank'.
    if (!is.numeric(splicings$exon_rank)
     || any(is.na(splicings$exon_rank)))
        stop("'splicings$exon_rank' must be an integer vector with no NAs")
    if (!is.integer(splicings$exon_rank))
        splicings$exon_rank <- as.integer(splicings$exon_rank)
    if (any(splicings$exon_rank <= 0L))
        stop("'splicings$exon_rank' contains non-positive values")
    ## Check 'exon_id'.
    if ("exon_id" %in% colnames(splicings)
     && (!is.integer(splicings$exon_id) || any(is.na(splicings$exon_id))))
        stop("'splicings$exon_id' must be an integer vector, with no NAs")
    ## Check 'exon_name'.
    if ("exon_name" %in% colnames(splicings)
     && !.isCharacterVectorOrFactor(splicings$exon_name))
        stop("'splicings$exon_name' must be a character vector (or factor)")
    ## Check 'exon_chrom'.
    if ("exon_chrom" %in% colnames(splicings)
     && (!.isCharacterVectorOrFactor(splicings$exon_chrom)
         || any(is.na(splicings$exon_chrom))))
        stop("'splicings$exon_chrom' must be a character vector (or factor) ",
             "with no NAs")
    ## Check 'exon_strand'.
    if ("exon_strand" %in% colnames(splicings)
     && (!.isCharacterVectorOrFactor(splicings$exon_strand)
         || any(is.na(splicings$exon_strand))))
        stop("'splicings$exon_strand' must be a character vector (or factor) ",
             "with no NAs")
    if (("exon_chrom" %in% colnames(splicings))
     && !("exon_strand" %in% colnames(splicings)))
        stop("if 'splicings' has an \"exon_chrom\" col then ",
             "it must have an \"exon_strand\" col too")
    ## Check 'exon_start'.
    if (!is.numeric(splicings$exon_start)
     || any(is.na(splicings$exon_start)))
        stop("'splicings$exon_start' must be an integer vector with no NAs")
    if (!is.integer(splicings$exon_start))
        splicings$exon_start <- as.integer(splicings$exon_start)
    ## Check 'exon_end'.
    if (!is.numeric(splicings$exon_end)
     || any(is.na(splicings$exon_end)))
        stop("'splicings$exon_end' must be an integer vector with no NAs")
    if (!is.integer(splicings$exon_end))
        splicings$exon_end <- as.integer(splicings$exon_end)
    ## Check 'exon_start <= exon_end'.
    if (any(splicings$exon_start > splicings$exon_end))
        stop("exon starts must be <= exon ends")
    ## Check 'cds_start', 'cds_end'.
    if (is.null(splicings$cds_start) != is.null(splicings$cds_end))
        stop("'splicings' has a \"cds_start\" col ",
             "but no \"cds_end\" col, or vice versa")
    if (is.null(splicings$cds_start)) {
        warning("no CDS information for this TranscriptDb object")
    } else {
        if (!is.numeric(splicings$cds_start))
            stop("'splicings$cds_start' must be an integer vector")
        if (!is.integer(splicings$cds_start))
            splicings$cds_start <- as.integer(splicings$cds_start)
        if (!is.numeric(splicings$cds_end))
            stop("'splicings$cds_end' must be an integer vector")
        if (!is.integer(splicings$cds_end))
            splicings$cds_end <- as.integer(splicings$cds_end)
        if (!all(is.na(splicings$cds_end) == is.na(splicings$cds_start)))
            stop("NAs in 'splicings$cds_end' don't match ",
                 "NAs in 'splicings$cds_start'")
        if (any(splicings$cds_start > splicings$cds_end, na.rm=TRUE))
            stop("cds starts must be <= cds ends")
        if (any(splicings$cds_start < splicings$exon_start, na.rm=TRUE)
         || any(splicings$cds_end > splicings$exon_end, na.rm=TRUE))
            stop("cds starts/ends are incompatible with exon starts/ends")
    }
    ## Check 'cds_id'.
    if (!is.null(splicings$cds_id)) {
        if (is.null(splicings$cds_start))
            stop("'splicings' has a \"cds_id\" col ",
                 "but no \"cds_start\"/\"cds_end\" cols")
        if (!is.integer(splicings$cds_id))
            stop("'splicings$cds_id' must be an integer vector")
        if (!all(is.na(splicings$cds_id) == is.na(splicings$cds_start)))
            stop("NAs in 'splicings$cds_id' don't match ",
                 "NAs in 'splicings$cds_start'")
    }
    ## Check 'cds_name'.
    if (!is.null(splicings$cds_name)) {
        if (is.null(splicings$cds_start))
            stop("'splicings' has a \"cds_name\" col ",
                 "but no \"cds_start\"/\"cds_end\" cols")
        if (!.isCharacterVectorOrFactor(splicings$cds_name))
            stop("'splicings$cds_name' must be a character vector (or factor)")
        if (!all(is.na(splicings$cds_name) == is.na(splicings$cds_start)))
            stop("NAs in 'splicings$cds_name' don't match ",
                 "NAs in 'splicings$cds_start'")
    }
    splicings
}

.normargGenes <- function(genes, unique_tx_ids)
{
    if (is.null(genes))
        return(data.frame(tx_id=unique_tx_ids[FALSE], gene_id=character(0)))
    .REQUIRED_COLS <- "gene_id"
    .OPTIONAL_COLS <- c("tx_id", "tx_name")
    .checkargColnames(genes, .REQUIRED_COLS, .OPTIONAL_COLS, "genes")
    ## Check 'gene_id'.
    if (!.isCharacterVectorOrFactor(genes$gene_id)
     || any(is.na(genes$gene_id)))
        stop("'genes$gene_id' must be a character vector (or factor) ",
             "with no NAs")
    ## 'genes' must have one of the 2 optional cols but not both.
    if (length(intersect(colnames(genes), .OPTIONAL_COLS)) != 1L)
        stop("'genes' must have either a \"tx_id\" ",
             "or a \"tx_name\" col but not both")
    if (is.null(genes$tx_id)) {
        ## Remap 'gene_id' to 'tx_id'.
        if (is.null(names(unique_tx_ids)))
            stop("cannot map genes to transcripts, need 'transcripts$tx_name'")
        genes <- joinDataFrameWithName2Val(genes, "tx_name",
                                           unique_tx_ids, "tx_id")
    } else {
        ## Check 'tx_id'.
        .checkForeignKey(genes$tx_id, "integer", "genes$tx_id",
                         unique_tx_ids, "integer", "transcripts$tx_id")
    }
    genes
}

.normargChrominfo <- function(chrominfo, transcripts_tx_chrom,
                              splicings_exon_chrom)
{
    if (is.null(chrominfo)) {
        feature_chrom <- unique(c(as.character(transcripts_tx_chrom),
                                  as.character(splicings_exon_chrom)))
        chrominfo <- data.frame(
            chrom=feature_chrom,
            length=rep.int(NA_integer_, length(feature_chrom))
        )
        return(chrominfo)
    }
    .REQUIRED_COLS <- c("chrom", "length")
    .OPTIONAL_COLS <- character(0)
    .checkargColnames(chrominfo, .REQUIRED_COLS, .OPTIONAL_COLS, "chrominfo")
    ## Check 'chrom'.
    if (!.isCharacterVectorOrFactor(chrominfo$chrom)
     || any(is.na(chrominfo$chrom)))
        stop("'chrominfo$chrom' must be a character vector (or factor) ",
             "with no NAs")
    .checkForeignKey(transcripts_tx_chrom, NA, "transcripts$tx_chrom",
                     chrominfo$chrom, NA, "chrominfo$chrom")
    if (!is.null(splicings_exon_chrom))
        .checkForeignKey(splicings_exon_chrom, NA, "splicings$exon_chrom",
                         chrominfo$chrom, NA, "chrominfo$chrom")
    ## Check 'length'.
    if (!is.vector(chrominfo$length))
        stop("'chrominfo$length' must be either all NAs ",
             "or an integer vector with no NAs")
    na_idx <- is.na(chrominfo$length)
    if (!all(na_idx)) {
        if (any(na_idx))
            stop("'chrominfo$length' cannot mix NAs and non-NAs")
        if (!is.numeric(chrominfo$length))
            stop("'chrominfo$length' must be either all NAs ",
                 "or an integer vector with no NAs")
    }
    if (!is.integer(chrominfo$length))
        chrominfo$length <- as.integer(chrominfo$length)
    chrominfo
}

.makeInternalIdsFromExternalIds <- function(external_id)
{
    if (is.integer(external_id))
        external_id
    else
        as.integer(factor(external_id))
}

.makeInternalIdsForUniqueLocs <- function(chrom, strand, start, end)
{
    not_NA <- !is.na(start)
    x <- data.frame(chrom, strand, start, end,
                    stringsAsFactors=FALSE)[not_NA, ]
    ans <- integer(length(start))
    ans[not_NA] <- makeIdsForUniqueDataFrameRows(x)
    ans[!not_NA] <- NA_integer_
    ans
}

.writeFeatureTable <- function(conn,
                               tablename,
                               internal_id,
                               name,
                               chrom,
                               strand,
                               start,
                               end,
                               feature_shortname=NA)
{
    if (is.null(name))
        name <- rep.int(NA_character_, length(internal_id))
    if (is.na(feature_shortname))
        feature_shortname <- tablename
    colnames <- .makeFeatureColnames(feature_shortname)
    table <- data.frame(
        internal_id=internal_id,
        name=name,
        chrom=chrom,
        strand=strand,
        start=start,
        end=end,
        stringsAsFactors=FALSE)
    table <- unique(table)

    ## Create the '<tablename>' table.
    sql <- c(
        "CREATE TABLE ", tablename, " (\n",
        "  ", colnames[1L], " INTEGER PRIMARY KEY,\n",
        "  ", colnames[2L], " TEXT NULL,\n",
        "  ", colnames[3L], " TEXT NOT NULL,\n",
        "  ", colnames[4L], " TEXT NOT NULL,\n",
        "  ", colnames[5L], " INTEGER NOT NULL,\n",
        "  ", colnames[6L], " INTEGER NOT NULL\n",
        ")")
    res <- dbSendQuery(conn, paste(sql, collapse=""))
    dbClearResult(res)

    ## Fill the '<tablename>' table.
    ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
    ## when the nb of rows to insert is 0, hence the following test.
    if (nrow(table) != 0L) {
        sql <- c("INSERT INTO ", tablename, " VALUES (?,?,?,?,?,?)")
        dbBeginTransaction(conn)
        res <- dbSendPreparedQuery(conn, paste(sql, collapse=""),
                                   table)
        dbClearResult(res)
        dbCommit(conn)
    }
}

.writeSplicingTable <- function(conn,
                                internal_tx_id,
                                exon_rank,
                                internal_exon_id,
                                internal_cds_id)
{
    table <- data.frame(
        internal_tx_id=internal_tx_id,
        exon_rank=exon_rank,
        internal_exon_id=internal_exon_id,
        internal_cds_id=internal_cds_id,
        stringsAsFactors=FALSE)
    table <- unique(table)

    ## Create the 'splicing' table and related indices.
    sql <- c(
        "CREATE TABLE splicing (\n",
        "  _tx_id INTEGER NOT NULL,\n",
        "  exon_rank INTEGER NOT NULL,\n",
        "  _exon_id INTEGER NOT NULL,\n",
        "  _cds_id INTEGER NULL,\n",
        "  UNIQUE (_tx_id, exon_rank),\n",
        "  FOREIGN KEY (_tx_id) REFERENCES transcript,\n",
        "  FOREIGN KEY (_exon_id) REFERENCES exon,\n",
        "  FOREIGN KEY (_cds_id) REFERENCES cds\n",
        ")")
    res <- dbSendQuery(conn, paste(sql, collapse=""))
    dbClearResult(res)
    sql <- c(
        "CREATE INDEX F_tx_id ON splicing (_tx_id);\n",
        "CREATE INDEX F_exon_id ON splicing (_exon_id);\n",
        "CREATE INDEX F_cds_id ON splicing (_cds_id)"
    )
    #Temporarily droped the indices.
    #res <- dbSendQuery(conn, paste(sql, collapse=""))
    #dbClearResult(res)

    ## Fill the 'splicing' table.
    ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
    ## when the nb of rows to insert is 0, hence the following test.
    if (nrow(table) != 0L) {
        sql <- "INSERT INTO splicing VALUES (?,?,?,?)"
        dbBeginTransaction(conn)
        res <- dbSendPreparedQuery(conn, sql, table)
        dbClearResult(res)
        dbCommit(conn)
    }
}

.writeGeneTable <- function(conn, gene_id, internal_tx_id)
{
    table <- data.frame(
        gene_id=gene_id,
        internal_tx_id=internal_tx_id,
        stringsAsFactors=FALSE)
    table <- unique(table)
    table <- table[!is.na(table$gene_id), ]
    ## Create the 'gene' table.
    sql <- c(
        "CREATE TABLE gene (\n",
        "  gene_id TEXT NOT NULL,\n",
        "  _tx_id INTEGER NOT NULL,\n",
        "  UNIQUE (gene_id, _tx_id),\n",
        "  FOREIGN KEY (_tx_id) REFERENCES transcript\n",
        ")")
    res <- dbSendQuery(conn, paste(sql, collapse=""))
    dbClearResult(res)
    ## Fill the 'gene' table.
    ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
    ## when the nb of rows to insert is 0, hence the following test.
    if (nrow(table) != 0L) {
        sql <- "INSERT INTO gene VALUES (?,?)"
        dbBeginTransaction(conn)
        res <- dbSendPreparedQuery(conn, sql, table)
        dbClearResult(res)
        dbCommit(conn)
    }
}

.writeChrominfoTable <- function(conn, chrominfo)
{
    table <- data.frame(
        internal_chrom_id=seq_len(nrow(chrominfo)),
        chrom=as.character(chrominfo$chrom),
        length=chrominfo$length,
        stringsAsFactors=FALSE)
    ## Create the 'chrominfo' table.
    sql <- c(
        "CREATE TABLE chrominfo (\n",
        "  _chrom_id INTEGER PRIMARY KEY,\n",
        "  chrom TEXT UNIQUE NOT NULL,\n",
        "  length INTEGER NULL\n",
        ")")
    res <- dbSendQuery(conn, paste(sql, collapse=""))
    dbClearResult(res)
    ## Fill the 'chrominfo' table.
    ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
    ## when the nb of rows to insert is 0, hence the following test.
    if (nrow(table) != 0L) {
        sql <- "INSERT INTO chrominfo VALUES (?,?,?)"
        dbBeginTransaction(conn)
        res <- dbSendPreparedQuery(conn, sql, table)
        dbClearResult(res)
        dbCommit(conn)
    }
}

.writeMetadataTable <- function(conn, metadata)
{
    transcript_nrow <- dbGetQuery(conn, "SELECT COUNT(*) FROM transcript")[[1L]]
    exon_nrow <- dbGetQuery(conn, "SELECT COUNT(*) FROM exon")[[1L]]
    cds_nrow <- dbGetQuery(conn, "SELECT COUNT(*) FROM cds")[[1L]]
    thispkg_version <- installed.packages()['GenomicFeatures', 'Version']
    rsqlite_version <- installed.packages()['RSQLite', 'Version']
    mat1 <- matrix(c(
        .DB_TYPE_NAME,     .DB_TYPE_VALUE),
        ncol=2, byrow=TRUE
    )
    mat2 <- matrix(c(
        "transcript_nrow", transcript_nrow,
        "exon_nrow",       exon_nrow,
        "cds_nrow",        cds_nrow,
        "Db created by",   "GenomicFeatures package from Bioconductor",
        "Creation time",   svn.time(),
        "GenomicFeatures version at creation time", thispkg_version,
        "RSQLite version at creation time", rsqlite_version),
        ncol=2, byrow=TRUE
    )
    colnames(mat1) <- colnames(mat2) <- c("name", "value")
    metadata <- rbind(data.frame(name=mat1[ , "name"], value=mat1[ , "value"],
                                 stringsAsFactors=FALSE),
                      metadata,
                      data.frame(name=mat2[ , "name"], value=mat2[ , "value"],
                                 stringsAsFactors=FALSE))
    dbWriteTable(conn, "metadata", metadata, row.names=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTranscriptDb().
###

.importTranscripts <- function(conn, transcripts, internal_tx_id)
{
    .writeFeatureTable(conn, "transcript",
        internal_tx_id, transcripts$tx_name,
        transcripts$tx_chrom, transcripts$tx_strand,
        transcripts$tx_start, transcripts$tx_end,
        feature_shortname="tx")
}

.importExons <- function(conn, splicings, internal_exon_id)
{
    .writeFeatureTable(conn, "exon",
        internal_exon_id, splicings$exon_name,
        splicings$exon_chrom, splicings$exon_strand,
        splicings$exon_start, splicings$exon_end)
}

.importCDS <- function(conn, splicings, internal_cds_id)
{
    cds_name <- splicings$cds_name[!is.na(internal_cds_id)]
    cds_chrom <- splicings$exon_chrom[!is.na(internal_cds_id)]
    cds_strand <- splicings$exon_strand[!is.na(internal_cds_id)]
    cds_start <- splicings$cds_start[!is.na(internal_cds_id)]
    cds_end <- splicings$cds_end[!is.na(internal_cds_id)]
    .writeFeatureTable(conn, "cds",
        internal_cds_id[!is.na(internal_cds_id)], cds_name,
        cds_chrom, cds_strand,
        cds_start, cds_end)
}

makeTranscriptDb <- function(transcripts, splicings,
                             genes=NULL, chrominfo=NULL, metadata=NULL, ...)
{
    if (length(list(...)) != 0L)
        warning("extra args are ignored for now")
    transcripts <- .normargTranscripts(transcripts)
    unique_tx_ids <- transcripts$tx_id  # guaranteed to be unique
    names(unique_tx_ids) <- transcripts$tx_name
    splicings <- .normargSplicings(splicings, unique_tx_ids)
    genes <- .normargGenes(genes, unique_tx_ids)
    chrominfo <- .normargChrominfo(chrominfo, transcripts$tx_chrom,
                                   splicings$exon_chrom)
    transcripts_internal_tx_id <- unique_tx_ids
    splicings_internal_tx_id <- splicings$tx_id
    genes_internal_tx_id <- genes$tx_id
    ## Infer 'splicings$exon_chrom' and 'splicings$exon_strand' when missing
    ## and generate internal exon id.
    splicings2transcripts <- match(splicings_internal_tx_id, unique_tx_ids)
    if (is.null(splicings$exon_chrom))
        splicings$exon_chrom <- transcripts$tx_chrom[splicings2transcripts]
    if (is.null(splicings$exon_strand))
        splicings$exon_strand <- transcripts$tx_strand[splicings2transcripts]
    if (!is.null(splicings$exon_id)) {
        splicings_internal_exon_id <- splicings$exon_id
    } else if (!is.null(splicings$exon_name)) {
        splicings_internal_exon_id <-
            .makeInternalIdsFromExternalIds(splicings$exon_name)
        #splicings$exon_id <- splicings_internal_exon_id
    } else {
        splicings_internal_exon_id <-
            .makeInternalIdsForUniqueLocs(
                splicings$exon_chrom, splicings$exon_strand,
                splicings$exon_start, splicings$exon_end)
        #splicings$exon_id <- splicings_internal_exon_id
    }
    ## Infer 'splicings$cds_start' and 'splicings$cds_end' when missing
    ## and generate internal cds id.
    if (is.null(splicings$cds_start)) {
        splicings$cds_start <- rep.int(NA_integer_, nrow(splicings))
        splicings$cds_end <- splicings$cds_start
    }
    if (!is.null(splicings$cds_id)) {
        splicings_internal_cds_id <- splicings$cds_id
    } else if (!is.null(splicings$cds_name)) {
        splicings_internal_cds_id <-
            .makeInternalIdsFromExternalIds(splicings$cds_name)
        #splicings$cds_id <- splicings_internal_cds_id
    } else {
        splicings_internal_cds_id <-
            .makeInternalIdsForUniqueLocs(
                splicings$exon_chrom, splicings$exon_strand,
                splicings$cds_start, splicings$cds_end)
        #splicings$cds_id <- splicings_internal_cds_id
    }
    ## Create the db in a temp file.
    conn <- dbConnect(SQLite(), dbname="")
    .importTranscripts(conn, transcripts, transcripts_internal_tx_id)
    .importExons(conn, splicings, splicings_internal_exon_id)
    .importCDS(conn, splicings, splicings_internal_cds_id)
    .writeSplicingTable(conn,
                        splicings_internal_tx_id,
                        splicings$exon_rank,
                        splicings_internal_exon_id,
                        splicings_internal_cds_id)
    .writeGeneTable(conn, genes$gene_id, genes_internal_tx_id)
    .writeChrominfoTable(conn, chrominfo)
    .writeMetadataTable(conn, metadata)  # must come last!
    new("TranscriptDb", conn=conn)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method for TranscriptDb objects.
###

setMethod("show", "TranscriptDb",
    function(object)
    {
        cat("TranscriptDb object:\n")
        metadata <- dbReadTable(object@conn, "metadata")
        for (i in seq_len(nrow(metadata))) {
            cat("| ", metadata[i, "name"], ": ", metadata[i, "value"],
                "\n", sep="")
        }
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Comparing 2 TranscriptDb objects.
###

### Dump the entire db into a list of data frames 'txdump' that can be used
### in 'do.call(makeTranscriptDb, txdump)' to make the db again with no loss
### of information.
### Note that the transcripts are dumped in the same order in all the
### data frames.
setMethod("as.list", "TranscriptDb",
    function(x, ...)
    {
        ORDER_BY <- "ORDER BY tx_chrom, tx_strand, tx_start, tx_end, tx_id"
        ## Retrieve the "transcripts" element.
        sql <- paste("SELECT transcript._tx_id AS tx_id, tx_name,",
                     "tx_chrom, tx_strand, tx_start, tx_end FROM transcript",
                     ORDER_BY)
        transcripts <- dbGetQuery(x@conn, sql)
        COL2CLASS <- c(
             tx_id="integer",
             tx_name="character",
             tx_chrom="factor",
             tx_strand="factor",
             tx_start="integer",
             tx_end="integer"
        )
        transcripts <- setDataFrameColClass(transcripts, COL2CLASS)

        ## Retrieve the "splicings" element.
        sql <- paste(
            "SELECT transcript._tx_id AS tx_id, exon_rank,",
            "exon._exon_id AS exon_id, exon_name,",
            "exon_chrom, exon_strand, exon_start, exon_end,",
            #"cds._cds_id AS cds_id, cds_name,",
            "cds._cds_id AS cds_id,",
            "cds_start, cds_end",
            "FROM transcript",
            "INNER JOIN splicing",
            "ON (transcript._tx_id=splicing._tx_id)",
            "INNER JOIN exon",
            "ON (splicing._exon_id=exon._exon_id)",
            "LEFT JOIN cds",
            "ON (splicing._cds_id=cds._cds_id)",
            ORDER_BY, ", exon_rank")
        splicings <- dbGetQuery(x@conn, sql)
        COL2CLASS <- c(
             tx_id="integer",
             exon_rank="integer",
             exon_id="integer",
             exon_name="character",
             exon_chrom="factor",
             exon_strand="factor",
             exon_start="integer",
             exon_end="integer",
             cds_id="integer",
             #cds_name="character",
             cds_start="integer",
             cds_end="integer"
        )
        splicings <- setDataFrameColClass(splicings, COL2CLASS)

        ## Retrieve the "genes" element.
        sql <- paste(
            "SELECT transcript._tx_id AS tx_id, gene_id",
            "FROM transcript",
            "INNER JOIN gene",
            "ON (transcript._tx_id=gene._tx_id)",
            ORDER_BY, ", gene_id")
        genes <- dbGetQuery(x@conn, sql)
        COL2CLASS <- c(
             tx_id="integer",
             gene_id="character"
        )
        genes <- setDataFrameColClass(genes, COL2CLASS)

        ## Retrieve the "chrominfo" element.
        sql <- "SELECT chrom, length FROM chrominfo ORDER BY _chrom_id"
        chrominfo <- dbGetQuery(x@conn, sql)

        list(transcripts=transcripts, splicings=splicings,
             genes=genes, chrominfo=chrominfo)
    }
)

compareTranscriptDbs <- function(txdb1, txdb2)
{
    if (!is(txdb1, "TranscriptDb")
     || !is(txdb2, "TranscriptDb"))
        stop("'txdb1' and 'txdb2' must be TranscriptDb objects")
    txdump1 <- as.list(txdb1)
    txdump2 <- as.list(txdb2)
    identical(txdump1, txdump2)
}

