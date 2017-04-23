### =========================================================================
### Making TxDb objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 1st group of helper functions for makeTxDb()
###
### 4 functions to check and normalize the input of makeTxDb():
###   o .normarg_makeTxDb_transcripts()
###   o .normarg_makeTxDb_splicings()
###   o .normarg_makeTxDb_genes()
###   o .normarg_makeTxDb_chrominfo()

.all_logical_NAs <- function(x)
{
    is.logical(x) && all(is.na(x))
}

### Also returns TRUE if 'x' is a logical vector of NAs.
.is_character_or_factor <- function(x)
{
    is.character(x) || is.factor(x) || .all_logical_NAs(x)
}

### Like as.integer(x) but fail if 'x' is not a numeric vector or a logical
### vector of NAs.
.graceful_as_integer <- function(x, x_what="x")
{
    if (!(is.numeric(x) || .all_logical_NAs(x)))
        stop(wmsg("'", x_what, "' must be either all NAs ",
                  "or an integer vector"))
    if (!is.integer(x))
        x <- as.integer(x)
    x
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

.normarg_makeTxDb_transcripts <- function(transcripts)
{
    .REQUIRED_COLS <- c("tx_id", "tx_chrom", "tx_strand", "tx_start", "tx_end")
    .OPTIONAL_COLS <- c("tx_name", "tx_type")
    check_colnames(transcripts, .REQUIRED_COLS, .OPTIONAL_COLS, "transcripts")
    ## Check 'tx_id'.
    if (!is.integer(transcripts$tx_id) || any(is.na(transcripts$tx_id)))
        stop("'transcripts$tx_id' must be an integer vector, with no NAs")
    if (any(duplicated(transcripts$tx_id)))
        stop("'transcripts$tx_id' contains duplicated values")
    ## Check 'tx_name'.
    if (has_col(transcripts, "tx_name")
     && !.is_character_or_factor(transcripts$tx_name))
        stop("'transcripts$tx_name' must be a character vector (or factor)")
    ## Check 'tx_type'.
    if (has_col(transcripts, "tx_type")
     && !.is_character_or_factor(transcripts$tx_type))
        stop("'transcripts$tx_type' must be a character vector (or factor)")
    ## Check 'tx_chrom'.
    if (!.is_character_or_factor(transcripts$tx_chrom)
     || any(is.na(transcripts$tx_chrom)))
        stop("'transcripts$tx_chrom' must be a character vector (or factor) ",
             "with no NAs")
    ## Check 'tx_strand'.
    if (!.is_character_or_factor(transcripts$tx_strand)
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

.normarg_makeTxDb_splicings <- function(splicings, transcripts_tx_id)
{
    .REQUIRED_COLS <- c("tx_id", "exon_rank", "exon_start", "exon_end")
    .OPTIONAL_COLS <- c("exon_id", "exon_name", "exon_chrom", "exon_strand",
                        "cds_id", "cds_name", "cds_start", "cds_end")
    check_colnames(splicings, .REQUIRED_COLS, .OPTIONAL_COLS, "splicings")
    ## Check 'tx_id'.
    .checkForeignKey(splicings$tx_id, "integer", "splicings$tx_id",
                     transcripts_tx_id, "integer", "transcripts$tx_id")
    ## Check 'exon_rank'.
    if (!is.numeric(splicings$exon_rank)
     || any(is.na(splicings$exon_rank)))
        stop("'splicings$exon_rank' must be an integer vector with no NAs")
    if (!is.integer(splicings$exon_rank))
        splicings$exon_rank <- as.integer(splicings$exon_rank)
    if (any(splicings$exon_rank <= 0L))
        stop("'splicings$exon_rank' contains non-positive values")
    ## Check 'exon_id'.
    if (has_col(splicings, "exon_id")
     && (!is.integer(splicings$exon_id) || any(is.na(splicings$exon_id))))
        stop("'splicings$exon_id' must be an integer vector, with no NAs")
    ## Check 'exon_name'.
    if (has_col(splicings, "exon_name")
     && !.is_character_or_factor(splicings$exon_name))
        stop("'splicings$exon_name' must be a character vector (or factor)")
    ## Check 'exon_chrom'.
    if (has_col(splicings, "exon_chrom")
     && (!.is_character_or_factor(splicings$exon_chrom)
         || any(is.na(splicings$exon_chrom))))
        stop("'splicings$exon_chrom' must be a character vector (or factor) ",
             "with no NAs")
    ## Check 'exon_strand'.
    if (has_col(splicings, "exon_strand")
     && (!.is_character_or_factor(splicings$exon_strand)
         || any(is.na(splicings$exon_strand))))
        stop("'splicings$exon_strand' must be a character vector (or factor) ",
             "with no NAs")
    if (has_col(splicings, "exon_chrom") && !has_col(splicings, "exon_strand"))
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
    ## Check presence of 'cds_start' and 'cds_end'.
    if (has_col(splicings, "cds_start") != has_col(splicings, "cds_end"))
        stop("'splicings' has a \"cds_start\" col ",
             "but no \"cds_end\" col, or vice versa")
    if (!has_col(splicings, "cds_start")) {
        warning("making a TxDb object without CDS information")
    } else {
        ## Check 'cds_start'.
        splicings$cds_start <- .graceful_as_integer(splicings$cds_start,
                                                    "splicings$cds_start")
        ## Check 'cds_end'.
        splicings$cds_end <- .graceful_as_integer(splicings$cds_end,
                                                  "splicings$cds_end")
        ## Check 'cds_start' and 'cds_end' compatibility.
        if (!all(is.na(splicings$cds_start) == is.na(splicings$cds_end)))
            stop("NAs in 'splicings$cds_start' don't match ",
                 "NAs in 'splicings$cds_end'")
        if (any(splicings$cds_start > splicings$cds_end, na.rm=TRUE))
            stop("cds starts must be <= cds ends")
        ## Check CDS and exon compatibility.
        if (any(splicings$cds_start < splicings$exon_start, na.rm=TRUE)
         || any(splicings$cds_end > splicings$exon_end, na.rm=TRUE))
            stop("cds starts/ends are incompatible with exon starts/ends")
    }
    ## Check 'cds_id'.
    if (has_col(splicings, "cds_id")) {
        if (!has_col(splicings, "cds_start"))
            stop("'splicings' has a \"cds_id\" col ",
                 "but no \"cds_start\"/\"cds_end\" cols")
        if (!is.integer(splicings$cds_id))
            stop("'splicings$cds_id' must be an integer vector")
        if (!all(is.na(splicings$cds_id) == is.na(splicings$cds_start)))
            stop("NAs in 'splicings$cds_id' don't match ",
                 "NAs in 'splicings$cds_start'")
    }
    ## Check 'cds_name'.
    if (has_col(splicings, "cds_name")) {
        if (!has_col(splicings, "cds_start"))
            stop("'splicings' has a \"cds_name\" col ",
                 "but no \"cds_start\"/\"cds_end\" cols")
        if (!.is_character_or_factor(splicings$cds_name))
            stop("'splicings$cds_name' must be a character vector (or factor)")
        if (any(is.na(splicings$cds_name) < is.na(splicings$cds_start)))
            stop("'splicings$cds_start' and 'splicings$cds_end' contain NAs ",
                 "where 'splicings$cds_name' doesn't")
    }
    if (!has_col(splicings, "cds_start")) {
        splicings$cds_start <- rep.int(NA_integer_, nrow(splicings))
        splicings$cds_end <- splicings$cds_start
    }
    splicings
}

.normarg_makeTxDb_genes <- function(genes, transcripts_tx_id)
{
    if (is.null(genes)) {
        genes <- data.frame(tx_id=transcripts_tx_id[FALSE],
                            gene_id=character(0),
                            check.names=FALSE, stringsAsFactors=FALSE)
        return(genes)
    }
    .REQUIRED_COLS <- "gene_id"
    .OPTIONAL_COLS <- c("tx_id", "tx_name")
    check_colnames(genes, .REQUIRED_COLS, .OPTIONAL_COLS, "genes")
    ## Check 'gene_id'.
    if (!.is_character_or_factor(genes$gene_id)
     || any(is.na(genes$gene_id)))
        stop("'genes$gene_id' must be a character vector (or factor) ",
             "with no NAs")
    ## 'genes' must have one of the 2 optional cols but not both.
    if (length(intersect(colnames(genes), .OPTIONAL_COLS)) != 1L)
        stop("'genes' must have either a \"tx_id\" ",
             "or a \"tx_name\" col but not both")
    if (!has_col(genes, "tx_id")) {
        ## Remap 'gene_id' to 'tx_id'.
        if (is.null(names(transcripts_tx_id)))
            stop("cannot map genes to transcripts, need 'transcripts$tx_name'")
        genes <- joinDataFrameWithName2Val(genes, "tx_name",
                                           transcripts_tx_id, "tx_id")
    } else {
        ## Check 'tx_id'.
        .checkForeignKey(genes$tx_id, "integer", "genes$tx_id",
                         transcripts_tx_id, "integer", "transcripts$tx_id")
    }
    genes
}

.normarg_makeTxDb_chrominfo <- function(chrominfo, transcripts_tx_chrom,
                                        splicings_exon_chrom)
{
    if (is.null(chrominfo)) {
        warning("chromosome lengths and circularity flags ",
                "are not available for this TxDb object")
        chrom <- unique(c(as.character(transcripts_tx_chrom),
                          as.character(splicings_exon_chrom)))
        chrom_ids <- rankSeqlevels(chrom)
        chrom[chrom_ids] <- chrom
        chrominfo <- data.frame(
            chrom=chrom,
            length=rep.int(NA_integer_, length(chrom)),
            is_circular=rep.int(NA, length(chrom)),
            check.names=FALSE, stringsAsFactors=FALSE
        )
        return(chrominfo)
    }
    .REQUIRED_COLS <- c("chrom", "length")
    .OPTIONAL_COLS <- "is_circular"
    check_colnames(chrominfo, .REQUIRED_COLS, .OPTIONAL_COLS, "chrominfo")
    ## Check 'chrom'.
    if (!.is_character_or_factor(chrominfo$chrom)
     || any(is.na(chrominfo$chrom)))
        stop("'chrominfo$chrom' must be a character vector (or factor) ",
             "with no NAs")
    if (any(duplicated(chrominfo$chrom)))
        stop("'chrominfo$chrom' contains duplicated values")
    .checkForeignKey(transcripts_tx_chrom, NA, "transcripts$tx_chrom",
                     chrominfo$chrom, NA, "chrominfo$chrom")
    if (!is.null(splicings_exon_chrom))
        .checkForeignKey(splicings_exon_chrom, NA, "splicings$exon_chrom",
                         chrominfo$chrom, NA, "chrominfo$chrom")
    ## Check 'length'.
    chrominfo$length <- .graceful_as_integer(chrominfo$length,
                                             "chrominfo$length")
    na_idx <- is.na(chrominfo$length)
    if (any(na_idx) && !all(na_idx))
        stop("'chrominfo$length' cannot mix NAs and non-NAs")
    ## Check 'is_circular'.
    if (!has_col(chrominfo, "is_circular")) {
        warning("chromosome circularity flags ",
                "are not available for this TxDb object")
        chrominfo$is_circular <- rep.int(NA, nrow(chrominfo))
    } else if (!is.logical(chrominfo$is_circular)) {
        stop("'chrominfo$is_circular' must be a logical vector")
    }
    chrominfo
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2nd group of helper functions for makeTxDb()
###
### These functions deal with id assignment and reassignment.
###

.makeTranscriptsInternalTxId <- function(transcripts, reassign.ids,
                                         chrominfo_chrom)
{
    if (!reassign.ids)
        return(transcripts$tx_id)
    chrom_ids <- match(transcripts$tx_chrom, chrominfo_chrom)
    makeFeatureIds(chrom_ids, transcripts$tx_strand,
                   transcripts$tx_start,
                   transcripts$tx_end,
                   name=transcripts$tx_name,
                   type=transcripts$tx_type)
}

.makeSplicingsInternalExonId <- function(splicings, reassign.ids,
                                         chrominfo_chrom)
{
    if (!reassign.ids && has_col(splicings, "exon_id"))
        return(splicings$exon_id)
    chrom_ids <- match(splicings$exon_chrom, chrominfo_chrom)
    makeFeatureIds(chrom_ids, splicings$exon_strand,
                   splicings$exon_start,
                   splicings$exon_end,
                   name=splicings$exon_name,
                   same.id.for.dups=TRUE)
}

.makeSplicingsInternalCDSId <- function(splicings, reassign.ids,
                                        chrominfo_chrom)
{
    if (!reassign.ids && has_col(splicings, "cds_id"))
        return(splicings$cds_id)
    chrom_ids <- match(splicings$exon_chrom, chrominfo_chrom)
    not_NA <- !is.na(splicings$cds_start)
    ids <- makeFeatureIds(chrom_ids[not_NA], splicings$exon_strand[not_NA],
                          splicings$cds_start[not_NA],
                          splicings$cds_end[not_NA],
                          name=splicings$cds_name[not_NA],
                          same.id.for.dups=TRUE)
    ans <- integer(nrow(splicings))
    ans[not_NA] <- ids
    ans[!not_NA] <- NA_integer_
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 3rd group of helper functions for makeTxDb()
###
### These functions deal with writing data to the database.
###

.writeChrominfoTable <- function(conn, chrominfo)
{
    table <- data.frame(
        internal_chrom_id=seq_len(nrow(chrominfo)),
        chrom=as.character(chrominfo$chrom),
        length=chrominfo$length,
        is_circular=chrominfo$is_circular,
        check.names=FALSE, stringsAsFactors=FALSE)
    ## Create the 'chrominfo' table.
    sql <- c(
        "CREATE TABLE chrominfo (\n",
        "  _chrom_id INTEGER PRIMARY KEY,\n",
        "  chrom TEXT UNIQUE NOT NULL,\n",
        "  length INTEGER NULL,\n",
        "  is_circular INTEGER NULL\n",
        ")")
    dbEasyQuery(conn, paste(sql, collapse=""))
    ## Fill the 'chrominfo' table.
    sql <- "INSERT INTO chrominfo VALUES (?,?,?,?)"
    dbEasyPreparedQuery(conn, sql, table)
}

.writeFeatureTable <- function(conn, tablename,
                               internal_id, name, type,
                               chrom, strand,
                               start, end,
                               feature_shortname=NA)
{
    if (is.na(feature_shortname))
        feature_shortname <- tablename
    colnames <- makeFeatureColnames(feature_shortname)
    if (is.null(name))
        name <- rep.int(NA_character_, length(internal_id))
    table <- data.frame(
        internal_id=internal_id,
        name=name,
        check.names=FALSE, stringsAsFactors=FALSE)
    if ("type" %in% names(colnames)) {
        if (is.null(type))
            type <- rep.int(NA_character_, length(internal_id))
        table$type <- type
    }
    table <- cbind(table,
                   data.frame(
                       chrom=chrom,
                       strand=strand,
                       start=start,
                       end=end,
                       check.names=FALSE, stringsAsFactors=FALSE))
    table <- unique(table)

    ## Create the '<tablename>' table.
    sql <- c(
        "CREATE TABLE ", tablename, " (\n",
        "  ", colnames[["id"]], " INTEGER PRIMARY KEY,\n",
        "  ", colnames[["name"]], " TEXT NULL,\n")
    if ("type" %in% names(colnames))
        sql <- c(sql, "  ", colnames[["type"]], " TEXT NULL,\n")
    sql <- c(sql,
        "  ", colnames[["chrom"]], " TEXT NOT NULL,\n",
        "  ", colnames[["strand"]], " TEXT NOT NULL,\n",
        "  ", colnames[["start"]], " INTEGER NOT NULL,\n",
        "  ", colnames[["end"]], " INTEGER NOT NULL,\n",
        "  FOREIGN KEY (", colnames[["chrom"]], ") REFERENCES chrominfo (chrom)\n",
        ")")
    dbEasyQuery(conn, paste(sql, collapse=""))

    ## Fill the '<tablename>' table.
    values <- paste0(rep.int("?", ncol(table)), collapse=",")
    sql <- c("INSERT INTO ", tablename, " VALUES (", values, ")")
    dbEasyPreparedQuery(conn, paste(sql, collapse=""), table)
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
        check.names=FALSE, stringsAsFactors=FALSE)
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
    dbEasyQuery(conn, paste(sql, collapse=""))
    sql <- c(
        "CREATE INDEX F_tx_id ON splicing (_tx_id);\n",
        "CREATE INDEX F_exon_id ON splicing (_exon_id);\n",
        "CREATE INDEX F_cds_id ON splicing (_cds_id)"
    )
    #Temporarily droped the indices.
    #dbEasyQuery(conn, paste(sql, collapse=""))

    ## Fill the 'splicing' table.
    sql <- "INSERT INTO splicing VALUES (?,?,?,?)"
    dbEasyPreparedQuery(conn, sql, table)
}

.writeGeneTable <- function(conn, gene_id, internal_tx_id)
{
    table <- data.frame(
        gene_id=gene_id,
        internal_tx_id=internal_tx_id,
        check.names=FALSE, stringsAsFactors=FALSE)
    table <- unique(table)
    table <- S4Vectors:::extract_data_frame_rows(table, !is.na(table$gene_id))
    ## Create the 'gene' table.
    sql <- c(
        "CREATE TABLE gene (\n",
        "  gene_id TEXT NOT NULL,\n",
        "  _tx_id INTEGER NOT NULL,\n",
        "  UNIQUE (gene_id, _tx_id),\n",
        "  FOREIGN KEY (_tx_id) REFERENCES transcript\n",
        ")")
    dbEasyQuery(conn, paste(sql, collapse=""))
    ## Fill the 'gene' table.
    sql <- "INSERT INTO gene VALUES (?,?)"
    dbEasyPreparedQuery(conn, sql, table)
}

.writeMetadataTable <- function(conn, metadata)
{
    transcript_nrow <- dbEasyQuery(conn,
                                   "SELECT COUNT(*) FROM transcript")[[1L]]
    exon_nrow <- dbEasyQuery(conn, "SELECT COUNT(*) FROM exon")[[1L]]
    cds_nrow <- dbEasyQuery(conn, "SELECT COUNT(*) FROM cds")[[1L]]
    thispkg_version <- packageDescription("GenomicFeatures")$Version
    rsqlite_version <- packageDescription("RSQLite")$Version
    mat1 <- matrix(c(
        DB_TYPE_NAME, DB_TYPE_VALUE,
        "Supporting package", "GenomicFeatures"),
        ncol=2, byrow=TRUE
    )
    mat2 <- matrix(c(
        "transcript_nrow", transcript_nrow,
        "exon_nrow",       exon_nrow,
        "cds_nrow",        cds_nrow,
        "Db created by",   "GenomicFeatures package from Bioconductor",
        "Creation time",   svn.time(),
        "GenomicFeatures version at creation time", thispkg_version,
        "RSQLite version at creation time", rsqlite_version,
        "DBSCHEMAVERSION", DB_SCHEMA_VERSION),
        ncol=2, byrow=TRUE
    )
    colnames(mat1) <- colnames(mat2) <- c("name", "value")
    metadata <- rbind(data.frame(name=mat1[ , "name"], value=mat1[ , "value"],
                                 check.names=FALSE, stringsAsFactors=FALSE),
                      metadata,
                      data.frame(name=mat2[ , "name"], value=mat2[ , "value"],
                                 check.names=FALSE, stringsAsFactors=FALSE))
    dbWriteTable(conn, "metadata", metadata, row.names=FALSE)
}

.importTranscripts <- function(conn, transcripts, internal_tx_id)
{
    .writeFeatureTable(conn, "transcript",
        internal_tx_id, transcripts$tx_name, transcripts$tx_type,
        transcripts$tx_chrom, transcripts$tx_strand,
        transcripts$tx_start, transcripts$tx_end,
        feature_shortname="tx")
}

.importExons <- function(conn, splicings, internal_exon_id)
{
    .writeFeatureTable(conn, "exon",
        internal_exon_id, splicings$exon_name, NULL,
        splicings$exon_chrom, splicings$exon_strand,
        splicings$exon_start, splicings$exon_end)
}

.importCDS <- function(conn, splicings, internal_cds_id)
{
    not_NA <- !is.na(internal_cds_id)
    internal_cds_id <- internal_cds_id[not_NA]
    cds_name <- splicings$cds_name[not_NA]
    cds_chrom <- splicings$exon_chrom[not_NA]
    cds_strand <- splicings$exon_strand[not_NA]
    cds_start <- splicings$cds_start[not_NA]
    cds_end <- splicings$cds_end[not_NA]
    .writeFeatureTable(conn, "cds",
        internal_cds_id, cds_name, NULL,
        cds_chrom, cds_strand,
        cds_start, cds_end)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTxDb().
###

makeTxDb <- function(transcripts, splicings,
                     genes=NULL, chrominfo=NULL, metadata=NULL,
                     reassign.ids=FALSE)
{
    if (!isTRUEorFALSE(reassign.ids))
        stop("'reassign.ids' must be TRUE or FALSE")
    transcripts <- .normarg_makeTxDb_transcripts(transcripts)
    transcripts_tx_id <- transcripts$tx_id  # guaranteed to be unique
    names(transcripts_tx_id) <- transcripts$tx_name
    splicings <- .normarg_makeTxDb_splicings(splicings, transcripts_tx_id)
    genes <- .normarg_makeTxDb_genes(genes, transcripts_tx_id)
    chrominfo <- .normarg_makeTxDb_chrominfo(chrominfo, transcripts$tx_chrom,
                                             splicings$exon_chrom)
    transcripts_internal_tx_id <- .makeTranscriptsInternalTxId(transcripts,
                                                               reassign.ids,
                                                               chrominfo$chrom)
    splicings_internal_tx_id <- translateIds(transcripts_tx_id,
                                             transcripts_internal_tx_id,
                                             splicings$tx_id)
    genes_internal_tx_id <- translateIds(transcripts_tx_id,
                                         transcripts_internal_tx_id,
                                         genes$tx_id)
    ## Infer 'splicings$exon_chrom' and 'splicings$exon_strand' when missing
    ## and generate internal exon id.
    splicings2transcripts <- match(splicings_internal_tx_id,
                                   transcripts_internal_tx_id)
    if (!has_col(splicings, "exon_chrom"))
        splicings$exon_chrom <- transcripts$tx_chrom[splicings2transcripts]
    if (!has_col(splicings, "exon_strand"))
        splicings$exon_strand <- transcripts$tx_strand[splicings2transcripts]
    splicings_internal_exon_id <- .makeSplicingsInternalExonId(splicings,
                                                               reassign.ids,
                                                               chrominfo$chrom)
    splicings_internal_cds_id <- .makeSplicingsInternalCDSId(splicings,
                                                             reassign.ids,
                                                             chrominfo$chrom)
    ## Create the db in a temp file.
    conn <- dbConnect(SQLite(), dbname="")
    .writeChrominfoTable(conn, chrominfo)  # must come first
    .importTranscripts(conn, transcripts, transcripts_internal_tx_id)
    .importExons(conn, splicings, splicings_internal_exon_id)
    .importCDS(conn, splicings, splicings_internal_cds_id)
    .writeSplicingTable(conn,
                        splicings_internal_tx_id,
                        splicings$exon_rank,
                        splicings_internal_exon_id,
                        splicings_internal_cds_id)
    .writeGeneTable(conn, genes$gene_id, genes_internal_tx_id)
    .writeMetadataTable(conn, metadata)  # must come last!
    TxDb(conn)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeToyTxDb().
###
### A simple wrapper around makeTxDb() typically used to make toy
### TxDb objects.
###
### 'splicings' must have the same format as for makeTxDb().
### If the "exon_chrom" col is missing, then it's added and filled with
### "chr1". If the "exon_strand" col is missing, then it's added and filled
### with "+". Within each transcript the "exon_chrom" and "exon_strand" values
### of the first and last exons must be the same (otherwise, there would be
### no easy way to infer the 'transcripts' table (see next).
###
### 'transcripts' is "inferred" from 'splicings' as follow:
###   - the "tx_chrom" and "tx_strand" cols are inferred from
###     'splicings$exon_chrom' and 'splicings$exon_strand' using the values
###     of the first exon (exon_rank=1) in the transcript;
###   - the transcripts starts/ends are inferred from the starts/ends of
###     their first and last exons.
###

makeToyTxDb <- function(splicings, genes=NULL)
{
    if (!is.data.frame(splicings))
        stop("'splicings' must be a data frame")
    stop("not ready yet, sorry!")
}


## helper to list mirbase.db miRBaseBuild values for species
supportedMiRBaseBuildValues <- function(){
    loadNamespace("mirbase.db")
    res <- toTable(mirbase.db::mirbaseSPECIES)[,c("name","genome_assembly")]
    S4Vectors:::extract_data_frame_rows(res, res$genome_assembly != "")
}
