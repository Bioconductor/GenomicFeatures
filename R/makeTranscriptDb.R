### =========================================================================
### Making TranscriptDb objects
### -------------------------------------------------------------------------


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

### FIXME: loadFeatures() needs to put the db back into memory (it's currently
### returning a TranscriptDb object that points to the on-disk db).
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

.checkIntForeignKey <- function(referring_vals, referred_vals,
                                referring_colname, referred_colname)
{
    if (!is.integer(referring_vals))
        stop("'", referring_colname, "' must be an integer vector")
    if (any(is.na(referring_vals)))
        stop("'", referring_colname, "' cannot contain NAs")
    if (!all(referring_vals %in% referred_vals))
        stop("all the values in '", referring_vals, "' must ",
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
    .checkIntForeignKey(splicings$tx_id, unique_tx_ids,
                        "splicings$tx_id", "transcripts$tx_id")
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
    if (!is.null(splicings$cds_start)) {
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
        .checkIntForeignKey(genes$tx_id, unique_tx_ids,
                            "genes$tx_id", "transcripts$tx_id")
    }
    genes
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

### Because we use SQLite "rtree" feature, .writeFeatureCoreTables() creates
### the 5 following tables:
###   (1) '<feature>'
###   (2) '<feature>_rtree'
###   (3) '<feature>_rtree_node'
###   (4) '<feature>_rtree_parent'
###   (5) '<feature>_rtree_rowid'
### Note that only (1) and (2) are explicitely created. (3), (4) and (5) are
### automatically created by the SQLite engine.
.writeFeatureCoreTables <- function(conn,
                                    feature,
                                    internal_id,
                                    name,
                                    chrom,
                                    strand,
                                    start,
                                    end,
                                    colnames)
{
    if (is.null(name))
        name <- rep.int(NA_character_, length(internal_id))
    table <- data.frame(
        internal_id=internal_id,
        name=name,
        chrom=chrom,
        strand=strand,
        start=start,
        end=end,
        stringsAsFactors=FALSE)
    table <- unique(table)

    ## Create the '<feature>' table.
    sql <- c(
        "CREATE TABLE ", feature, " (\n",
        "  ", colnames[1L], " INTEGER PRIMARY KEY,\n",
        "  ", colnames[2L], " TEXT NULL,\n",
        "  ", colnames[3L], " TEXT NOT NULL,\n",
        "  ", colnames[4L], " TEXT NOT NULL\n",
        ")")
    res <- dbSendQuery(conn, paste(sql, collapse=""))
    dbClearResult(res)

    ## Fill the '<feature>' table.
    ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
    ## when the nb of rows to insert is 0, hence the following test.
    if (nrow(table) != 0L) {
        sql <- c("INSERT INTO ", feature, " VALUES (?,?,?,?)")
        dbBeginTransaction(conn)
        res <- dbSendPreparedQuery(conn, paste(sql, collapse=""),
                                   table[1:4])
        dbClearResult(res)
        dbCommit(conn)
    }

    ## Create the '<feature>_rtree' table.
    sql <- c(
        "CREATE VIRTUAL TABLE ", feature, "_rtree USING rtree (\n",
        "  ", colnames[1L], " INTEGER PRIMARY KEY,\n",
        "  ", colnames[5L], " INTEGER, -- NOT NULL is implicit in rtree\n",
        "  ", colnames[6L], " INTEGER  -- NOT NULL is implicit in rtree\n",
        "  -- FOREIGN KEY (", colnames[1L], ") REFERENCES ", feature, "\n",
        ")")
    res <- dbSendQuery(conn, paste(sql, collapse=""))
    dbClearResult(res)

    ## Fill the '<feature>_rtree' table.
    if (nrow(table) != 0L) {
        sql <- c("INSERT INTO ", feature, "_rtree VALUES (?,?,?)")
        dbBeginTransaction(conn)
        res <- dbSendPreparedQuery(conn, paste(sql, collapse=""),
                                   table[c(1L, 5:6)])
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

.writeMetadataTable <- function(conn)
{
    transcript_nrow <- dbGetQuery(conn, "SELECT COUNT(*) FROM transcript")[[1L]]
    exon_nrow <- dbGetQuery(conn, "SELECT COUNT(*) FROM exon")[[1L]]
    cds_nrow <- dbGetQuery(conn, "SELECT COUNT(*) FROM cds")[[1L]]
    mat <- matrix(c(
        "DbType",  "TranscriptDb",
        "CreationDate", date(),
        "transcript_nrow", transcript_nrow,
        "exon_nrow", exon_nrow,
        "cds_nrow", cds_nrow),
        ncol=2, byrow=TRUE
    )
    colnames(mat) <- c("name", "value")
    metadata <- data.frame(name=mat[ , "name"], value=mat[ , "value"],
                           stringsAsFactors=FALSE)
    dbWriteTable(conn, "metadata", metadata, row.names=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTranscriptDb().
###

.importTranscripts <- function(conn, transcripts, internal_tx_id)
{
    .writeFeatureCoreTables(conn, "transcript",
        internal_tx_id, transcripts$tx_name,
        transcripts$tx_chrom, transcripts$tx_strand,
        transcripts$tx_start, transcripts$tx_end,
        c("_tx_id", "tx_name",
          "tx_chrom", "tx_strand",
          "tx_start", "tx_end"))
}

.importExons <- function(conn, splicings, internal_exon_id)
{
    .writeFeatureCoreTables(conn, "exon",
        internal_exon_id, splicings$exon_name,
        splicings$exon_chrom, splicings$exon_strand,
        splicings$exon_start, splicings$exon_end,
        c("_exon_id", "exon_name",
          "exon_chrom", "exon_strand",
          "exon_start", "exon_end"))
}

.importCDS <- function(conn, splicings, internal_cds_id)
{
    cds_name <- splicings$cds_name[!is.na(internal_cds_id)]
    cds_chrom <- splicings$exon_chrom[!is.na(internal_cds_id)]
    cds_strand <- splicings$exon_strand[!is.na(internal_cds_id)]
    cds_start <- splicings$cds_start[!is.na(internal_cds_id)]
    cds_end <- splicings$cds_end[!is.na(internal_cds_id)]
    .writeFeatureCoreTables(conn, "cds",
        internal_cds_id[!is.na(internal_cds_id)], cds_name,
        cds_chrom, cds_strand,
        cds_start, cds_end,
        c("_cds_id", "cds_name",
          "cds_chrom", "cds_strand",
          "cds_start", "cds_end"))
}

makeTranscriptDb <- function(transcripts, splicings, genes=NULL, ...)
{
    if (length(list(...)) != 0L)
        warning("extra args are ignored for now")
    transcripts <- .normargTranscripts(transcripts)
    unique_tx_ids <- transcripts$tx_id  # guaranteed to be unique
    names(unique_tx_ids) <- transcripts$tx_name
    splicings <- .normargSplicings(splicings, unique_tx_ids)
    genes <- .normargGenes(genes, unique_tx_ids)
    transcripts_internal_tx_id <- unique_tx_ids
    splicings_internal_tx_id <- splicings$tx_id
    genes_internal_tx_id <- genes$tx_id
    ## Infer 'splicings$exon_chrom' and 'splicings$exon_strand' when missing
    ## and generate internal exon id.
    if (is.null(splicings$exon_chrom))
        splicings$exon_chrom <- transcripts$tx_chrom[splicings_internal_tx_id]
    if (is.null(splicings$exon_strand))
        splicings$exon_strand <- transcripts$tx_strand[splicings_internal_tx_id]
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
    .writeMetadataTable(conn)
    new("TranscriptDb", conn=conn)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTranscriptDbFromUCSC() and UCSC-specific helper functions.
###

.makeExonRank <- function(exonCount, exonStrand)
{
    ans <- lapply(seq_len(length(exonCount)),
        function(i)
        {
            if (exonStrand[i] == "+")
                seq_len(exonCount[i])
            else
                (exonCount[i]):1L
        }
    )
    unlist(ans)
}

.extractExonRangesFromUCSCTxTable <- function(ucsc_txtable)
{
    exon_count <- ucsc_txtable$exonCount
    if (length(exon_count) == 0L)
        return(IRanges())
    ucsc_exonStarts <- strsplitAsListOfIntegerVectors(ucsc_txtable$exonStarts)
    if (!identical(elementLengths(ucsc_exonStarts), exon_count))
        stop("UCSC data anomaly: 'ucsc_txtable$exonStarts' ",
             "inconsistent with 'ucsc_txtable$exonCount'")
    exon_start <- unlist(ucsc_exonStarts) + 1L
    ucsc_exonEnds <- strsplitAsListOfIntegerVectors(ucsc_txtable$exonEnds)
    if (!identical(elementLengths(ucsc_exonEnds), exon_count))
        stop("UCSC data anomaly: 'ucsc_txtable$exonEnds' ",
             "inconsistent with 'ucsc_txtable$exonCount'")
    exon_end <- unlist(ucsc_exonEnds)
    IRanges(start=exon_start, end=exon_end)
}

### 'cds_start', 'cds_end': single integers
### 'exons': Ranges object
### Returns a list with 2 elements (start and end), each of them being integer
### vectors of the same length as 'exons', and may contain NAs.
### Notes:
###   (1) In refGene table, transcript NM_001146685: cds cumulative length is
###       not a multiple of 3:
###                 name chrom strand txStart   txEnd cdsStart  cdsEnd
###         NM_001146685  chr1      + 1351370 1353029  1351370 1353029
###         exonCount       exonStarts         exonEnds id   name2
###                 2 1351370,1352796, 1351628,1353029,  0 TMEM88B     
###         cdsStartStat cdsEndStat exonFrames
###                 cmpl     incmpl       0,0,
###       --> cds lengths: 1351628 - 1351370 -> 258
###                        1353029 - 1352796 -> 233
###       --> cds cum length: 491
###       Note that the cds end is marked as "incomplete" (see the cdsEndStat
###       col) which, according to UCSC, means that "the CDS is NOT completely
###       contained in the alignment at this end". See this post on the Genome
###       mailing list for more information:
###       https://lists.soe.ucsc.edu/pipermail/genome/2005-December/009184.html
###       Note that the post is about the Gencode Genes. Is it reasonable to
###       assume that this applies to RefSeq Genes too?
###   (2) Same thing in ensGene table, transcript ENST00000371841.
###   (3) All transcripts in knowGene table have a cds cumulative length that
###       is a multiple of 3.
### TODO: Investigate (1) and (2).
.extractUCSCCdsStartEnd <- function(cds_start, cds_end, exons, tx_name)
{
    ans_start <- ans_end <- integer(length(exons))
    ans_start[] <- NA_integer_
    ans_end[] <- NA_integer_
    if (cds_start > cds_end)
        return(list(start=ans_start, end=ans_end))
    first_exon_with_cds <- which(start(exons) <= cds_start
                                 & cds_start <= end(exons))
    if (length(first_exon_with_cds) != 1L)
        stop("UCSC data ambiguity in transcript ", tx_name,
             ": cannot determine first exon with cds ('cdsStart' ",
             "falls in 0 or more than 1 exon)")
    last_exon_with_cds <- which(start(exons) <= cds_end
                                & cds_end <= end(exons))
    if (length(last_exon_with_cds) != 1L)
        stop("UCSC data ambiguity in transcript ", tx_name,
             ": cannot determine last exon with cds ('cdsEnd' ",
             "falls in 0 or more than 1 exon)")
    if (last_exon_with_cds < first_exon_with_cds)
        stop("UCSC data anomaly in transcript ", tx_name,
             ": last exon with cds occurs before first exon with cds")
    exons_with_cds <- first_exon_with_cds:last_exon_with_cds
    ans_start[exons_with_cds] <- start(exons)[exons_with_cds]
    ans_end[exons_with_cds] <- end(exons)[exons_with_cds]
    ans_start[first_exon_with_cds] <- cds_start
    ans_end[last_exon_with_cds] <- cds_end
    if (sum(ans_end - ans_start + 1L, na.rm=TRUE) %% 3L != 0L)
        warning("UCSC data anomaly in transcript ", tx_name,
                ": the cds cumulative length is not a multiple of 3")
    list(start=ans_start, end=ans_end)
}

### Takes more than 5 minutes on the knownGene table! FIX IT!
.extractCdsRangesFromUCSCTxTable <- function(ucsc_txtable, exon_ranges)
{
    exon_count <- ucsc_txtable$exonCount
    cdsStart <- ucsc_txtable$cdsStart + 1L
    cdsEnd <- ucsc_txtable$cdsEnd
    cdsStarts <- cdsEnds <- vector(mode="list", length=nrow(ucsc_txtable))
    exon_ranges_offset <- 0L
    for (i in seq_len(nrow(ucsc_txtable))) {
        exons <- exon_ranges[exon_ranges_offset + seq_len(exon_count[i])]
        exon_ranges_offset <- exon_ranges_offset + exon_count[i]
        startend <- .extractUCSCCdsStartEnd(cdsStart[i], cdsEnd[i],
                                            exons, ucsc_txtable$name[i])
        cdsStarts[[i]] <- startend$start
        cdsEnds[[i]] <- startend$end
    }
    if (nrow(ucsc_txtable) == 0L) {
        ans_start <- integer(0)
        ans_end <- integer(0)
    } else {
        ans_start <- unlist(cdsStarts)
        ans_end <- unlist(cdsEnds)
    }
    return(list(start=ans_start, end=ans_end))
}

.makeTranscriptDbFromUCSCTxTable <- function(ucsc_txtable, genes, gene_id_type)
{
    COL2CLASS <- c(
        name="character",
        chrom="factor",
        strand="factor",
        txStart="integer",
        txEnd="integer",
        cdsStart="integer",
        cdsEnd="integer",
        exonCount="integer",
        exonStarts="character",
        exonEnds="character"
    )
    ucsc_txtable <- setDataFrameColClass(ucsc_txtable, COL2CLASS,
                                         drop.extra.cols=TRUE)

    ## Prepare the 'transcripts' data frame.
    tx_id <- seq_len(nrow(ucsc_txtable))
    transcripts <- data.frame(
        tx_id=tx_id,
        tx_name=ucsc_txtable$name,
        tx_chrom=ucsc_txtable$chrom,
        tx_strand=ucsc_txtable$strand,
        tx_start=ucsc_txtable$txStart + 1L,
        tx_end=ucsc_txtable$txEnd
    )

    ## Prepare the 'splicings' data frame.
    ## Exon starts and ends are multi-valued fields (comma-separated) that
    ## need to be expanded.
    exon_count <- ucsc_txtable$exonCount
    if (min(exon_count) <= 0L)
        stop("UCSC data anomaly: 'ucsc_txtable$exonCount' contains ",
             "non-positive values")
    splicings_tx_id <- rep.int(tx_id, exon_count)
    exon_ranges <- .extractExonRangesFromUCSCTxTable(ucsc_txtable)
    cds_startend <- .extractCdsRangesFromUCSCTxTable(ucsc_txtable, exon_ranges)
    splicings <- data.frame(
        tx_id=splicings_tx_id,
        exon_rank=.makeExonRank(exon_count, ucsc_txtable$strand),
        exon_start=start(exon_ranges),
        exon_end=end(exon_ranges),
        cds_start=cds_startend$start,
        cds_end=cds_startend$end
    )

    ## Prepare the 'genes' data frame.
    #genes <- genes[genes$tx_name %in% ucsc_txtable$name, ]

    ## Call makeTranscriptDb().
    makeTranscriptDb(transcripts, splicings, genes)
}

### Lookup between UCSC tables and tracks in the "Genes and Gene Prediction"
### group that are compatible with makeTranscriptDbFromUCSC(). A table is
### compatible if it has the following cols: name, chrom, strand, txStart,
### txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds.
### Note that from a strictly technical point of view the name and exonCount
### cols are not required (i.e. .makeTranscriptDbFromUCSCTxTable() could
### easily be modified to work even when they are missing).
.SUPPORTED_UCSC_TABLES <- c(
  ## tablename (unique key)    track             subtrack
  "knownGene",                 "UCSC Genes",     NA,
  "knownGeneOld3",             "Old UCSC Genes", NA,
  "wgEncodeGencodeManualRel2", "Gencode Genes",  "Genecode Manual",
  "wgEncodeGencodeAutoRel2",   "Gencode Genes",  "Genecode Auto",
  "wgEncodeGencodePolyaRel2",  "Gencode Genes",  "Genecode PolyA",
  "ccdsGene",                  "Consensus CDS",  NA, 
  "refGene",                   "RefSeq Genes",   NA, 
  "xenoRefGene",               "Other RefSeq",   NA, 
  "vegaGene",                  "Vega Genes",     "Vega Protein Genes", 
  "vegaPseudoGene",            "Vega Genes",     "Vega Pseudogenes", 
  "ensGene",                   "Ensembl Genes",  NA, 
  "acembly",                   "AceView Genes",  NA, 
  "sibGene",                   "SIB Genes",      NA, 
  "nscanPasaGene",             "N-SCAN",         "N-SCAN PASA-EST",
  "nscanGene",                 "N-SCAN",         "N-SCAN", 
  "sgpGene",                   "SGP Genes",      NA,
  "geneid",                    "Geneid Genes",   NA, 
  "genscan",                   "Genscan Genes",  NA, 
  "exoniphy",                  "Exoniphy",       NA, 
  "augustusHints",             "Augustus",       "Augustus Hints", 
  "augustusXRA",               "Augustus",       "Augustus De Novo", 
  "augustusAbinitio",          "Augustus",       "Augustus Ab Initio",
  "acescan",                   "ACEScan",        NA
)

supportedUCSCtables <- function()
{
    mat <- matrix(.SUPPORTED_UCSC_TABLES, ncol=3, byrow=TRUE)
    colnames(mat) <- c("tablename", "track", "subtrack")
    data.frame(track=mat[ , "track"], subtrack=mat[ , "subtrack"],
               row.names=mat[ , "tablename"],
               stringsAsFactors=FALSE)
}

### The table names above (unique key) must be used to name the elements of
### the list below. When a table name is missing, it means that the
### tx_name-to-gene_id mapping is not available so makeTranscriptDbFromUCSC()
### will leave the gene table empty.
.UCSC_TXNAME2GENEID_MAPINFO <- list(
    knownGene=c(             # left table (i.e. on the left side of the join)
        "knownToLocusLink",  # right table (i.e. on the right side of the join)
        "name",              # joining col in the right table
        "value",             # col in the right table that contains the gene id
        "Entrez Gene ID"),   # type of gene id
    wgEncodeGencodeManualRel2=c(
        "wgEncodeGencodeClassesRel2",
        "name",
        "geneId",
        "HAVANA Pseudogene ID"),
    wgEncodeGencodeAutoRel2=c(
        "wgEncodeGencodeClassesRel2",
        "name",
        "geneId",
        "HAVANA Pseudogene ID"),
    wgEncodeGencodePolyaRel2=c(
        "wgEncodeGencodeClassesRel2",
        "name",
        "geneId",
        "HAVANA Pseudogene ID"),
    refGene=c(
        "refLink",
        "mrnaAcc",
        "locusLinkId",
        "Entrez Gene ID"),
    vegaGene=c(
        "vegaGtp",
        "transcript",
        "gene",
        "HAVANA Pseudogene ID"),
    vegaPseudoGene=c(
        "vegaGtp",
        "transcript",
        "gene",
        "HAVANA Pseudogene ID"),
    ensGene=c(
        "ensGtp",
        "transcript",
        "gene",
        "Ensembl gene ID")
)

### The 2 main tasks that makeTranscriptDbFromUCSC() performs are:
###   (1) download the data from UCSC into a data.frame (the getTable() call);
###   (2) store that data.frame in an SQLite db (the
###       .makeTranscriptDbFromUCSCTxTable() call).
### Speed:
###   - for genome="hg18" and tablename="knownGene":
###       (1) download takes about 40-50 sec.
###       (2) db creation takes about 30-35 sec.
makeTranscriptDbFromUCSC <- function(genome="hg18", tablename="knownGene")
{
    if (!isSingleString(genome))
        stop("'genome' must be a single string")
    if (!isSingleString(tablename))
        stop("'tablename' must be a single string")
    track <- supportedUCSCtables()[tablename, "track"]
    if (is.na(track))
        stop("track \"", track, "\" is not supported")
    session <- browserSession()
    genome(session) <- genome
    ## Download the transcript table.
    query1 <- ucscTableQuery(session, track, table=tablename)
    ucsc_txtable <- getTable(query1)
    ## Download the tx_name-to-gene_id mapping.
    txname2gene_mapinfo <- .UCSC_TXNAME2GENEID_MAPINFO[[tablename]]
    if (is.null(txname2gene_mapinfo)) {
        genes <- NULL
        gene_id_type <- NA
    } else {
        tablename2 <- txname2gene_mapinfo[1L]
        query2 <- ucscTableQuery(session, track, table=tablename2)
        ucsc_genetable <- getTable(query2)
        tx_name <- ucsc_genetable[[txname2gene_mapinfo[2L]]]
        gene_id <- ucsc_genetable[[txname2gene_mapinfo[3L]]]
        if (is.null(tx_name) || is.null(gene_id))
            stop("expected cols \"", txname2gene_mapinfo[2L], "\" or/and \"",
                 txname2gene_mapinfo[3L], "\" not found in table ", tablename2)
        if (!is.character(tx_name))
            tx_name <- as.character(tx_name)
        if (!is.character(gene_id))
            gene_id <- as.character(gene_id)
        genes <- data.frame(tx_name=tx_name, gene_id=gene_id)
        gene_id_type <- txname2gene_mapinfo[4L]
    }
    .makeTranscriptDbFromUCSCTxTable(ucsc_txtable, genes, gene_id_type)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTranscriptDbFromBiomart().
###
### For people who want to tap BioMart.
### Typical use:
###   txdb <- makeTranscriptDbFromBiomart(biomart="ensembl",
###                                       dataset="hsapiens_gene_ensembl")
### Speed:
###   - for biomart="ensembl" and dataset="hsapiens_gene_ensembl":
###       (1) download takes about 8 min.
###       (2) db creation takes about 60-65 sec.
###

.extractCdsRangesFromBiomartTable <- function(bm_table)
{
    strand <- bm_table[["strand"]]
    cds_start <- exon_start <- bm_table[["exon_chrom_start"]]
    cds_end <- exon_end <- bm_table[["exon_chrom_end"]]
    utr5_start <- bm_table[["5_utr_start"]]
    utr5_end <- bm_table[["5_utr_end"]]
    utr3_start <- bm_table[["3_utr_start"]]
    utr3_end <- bm_table[["3_utr_end"]]

    if (!all(strand %in% c(1, -1)))
        stop("BioMart data anomaly: \"strand\" attribute should be 1 or -1")
    if (!is.numeric(exon_start)
     || !is.numeric(exon_end)
     || !is.numeric(utr5_start)
     || !is.numeric(utr5_end)
     || !is.numeric(utr3_start)
     || !is.numeric(utr3_end))
        stop("BioMart data anomaly: exon or utr coordinates don't ",
             "have a numeric type")
    no_utr5 <- is.na(utr5_start)
    if (!identical(no_utr5, is.na(utr5_end)))
        stop("BioMart data anomaly: NAs in \"5_utr_start\" attribute ",
             "don't match NAs in \"5_utr_end\" attribute")
    if (!all(utr5_start <= utr5_end, na.rm=TRUE))
        stop("BioMart data anomaly: some 5' UTR have a start > end")
    if (!all(utr5_start >= exon_start, na.rm=TRUE)
     || !all(utr5_end <= exon_end, na.rm=TRUE))
        stop("BioMart data anomaly: some 5' UTR are not within the exon limits")
    no_utr3 <- is.na(utr3_start)
    if (!identical(no_utr3, is.na(utr3_end)))
        stop("BioMart data anomaly: NAs in \"3_utr_start\" attribute ",
             "don't match NAs in \"3_utr_end\" attribute")
    if (!all(utr3_start <= utr3_end, na.rm=TRUE))
        stop("BioMart data anomaly: some 3' UTR have a start > end")
    if (!all(utr3_start >= exon_start, na.rm=TRUE)
     || !all(utr3_end <= exon_end, na.rm=TRUE))
        stop("BioMart data anomaly: some 3' UTR are not within the exon limits")

    idx <- strand == 1 & !no_utr5
    if (!all(utr5_start[idx] == exon_start[idx]))
        stop("BioMart data anomaly: some 5' UTR on the plus strand ",
             "don't start where the exon starts")
    cds_start[idx] <- utr5_end[idx] + 1L
    idx <- strand == 1 & !no_utr3
    if (!all(utr3_end[idx] == exon_end[idx]))
        stop("BioMart data anomaly: some 3' UTR on the plus strand ",
             "don't end where the exon ends")
    cds_end[idx] <- utr3_start[idx] - 1L
    idx <- strand == -1 & !no_utr3
    if (!all(utr3_start[idx] == exon_start[idx]))
        stop("BioMart data anomaly: some 3' UTR on the minus strand ",
             "don't start where the exon starts")
    cds_start[idx] <- utr3_end[idx] + 1L
    idx <- strand == -1 & !no_utr5
    if (!all(utr5_end[idx] == exon_end[idx]))
        stop("BioMart data anomaly: some 5' UTR on the minus strand ",
             "don't end where the exon ends")
    cds_end[idx] <- utr5_start[idx] - 1L
    ans <- IRanges(start=cds_start, end=cds_end)
    if (length(ans) != 0L) {
        cds_cumlength <-
            sapply(split(width(ans), bm_table$ensembl_transcript_id), sum)
        if (!all(cds_cumlength[as.vector(bm_table$ensembl_transcript_id)]
                 == bm_table$cds_length, na.rm=TRUE))
            stop("BioMart data anomaly: for some transcripts, the cds ",
                 "cumulative length inferred from the exon and UTR info ",
                 "doesn't match the \"cds_length\" attribute from BioMart")
        #if (!all(cds_cumlength %% 3L == 0L))
        #    warning("BioMart data anomaly: for some transcripts, the cds ",
        #            "cumulative length (\"cds_length\" attribute) is not ",
        #            "a multiple of 3")
    }
    ans
}

makeTranscriptDbFromBiomart <- function(biomart="ensembl",
                                        dataset="hsapiens_gene_ensembl",
                                        ensembl_transcript_ids=NULL)
{
    mart <- useMart(biomart=biomart, dataset=dataset)
    if (is.null(ensembl_transcript_ids)) {
        filters <- values <- ""
    } else if (is.character(ensembl_transcript_ids)
            && !any(is.na(ensembl_transcript_ids))) {
        filters <- "ensembl_transcript_id"
        values <- ensembl_transcript_ids
    } else {
            stop("'ensembl_transcript_ids' must be ",
                 "a character vector with no NAs")
    }
    ## Download and prepare the 'transcripts' data frame.
    attributes <- c("ensembl_transcript_id",
                    "chromosome_name",
                    "strand",
                    "transcript_start",
                    "transcript_end")
    bm_table <- getBM(attributes, filters=filters, values=values, mart=mart)
    transcripts_tx_id <- seq_len(nrow(bm_table))
    transcripts_tx_name <- bm_table$ensembl_transcript_id
    if (any(duplicated(transcripts_tx_name)))
        stop("the 'ensembl_transcript_id' field contains duplicates")
    transcripts <- data.frame(
        tx_id=transcripts_tx_id,
        tx_name=transcripts_tx_name,
        tx_chrom=bm_table$chromosome_name,
        tx_strand=ifelse(bm_table$strand == 1, "+", "-"),
        tx_start=bm_table$transcript_start,
        tx_end=bm_table$transcript_end
    )
    ## Download and prepare the 'splicings' data frame.
    ## Ironically the cds_start and cds_end attributes that we get from
    ## BioMart are pretty useless because they are relative to the coding
    ## mRNA. However, the utr coordinates are relative to the chromosome so
    ## we use them to infer the cds coordinates. We also retrieve the
    ## cds_length attribute as a sanity check.
    attributes <- c(
        "ensembl_transcript_id",
        "strand",
        "rank",
        "ensembl_exon_id",
        "exon_chrom_start",
        "exon_chrom_end",
        "5_utr_start",
        "5_utr_end",
        "3_utr_start",
        "3_utr_end",
        #"cds_start",
        #"cds_end",
        "cds_length"
    )
    bm_table <- getBM(attributes, filters=filters, values=values, mart=mart)
    splicings_tx_id <- as.integer(factor(bm_table$ensembl_transcript_id,
                                         levels=transcripts_tx_name))
    cds_ranges <- .extractCdsRangesFromBiomartTable(bm_table)
    cds_start <- start(cds_ranges)
    cds_start[width(cds_ranges) == 0L] <- NA_integer_
    cds_end <- end(cds_ranges)
    cds_end[width(cds_ranges) == 0L] <- NA_integer_
    splicings <- data.frame(
        tx_id=splicings_tx_id,
        exon_rank=bm_table$rank,
        exon_name=bm_table$ensembl_exon_id,
        exon_start=bm_table$exon_chrom_start,
        exon_end=bm_table$exon_chrom_end,
        cds_start=cds_start,
        cds_end=cds_end
    )
    ## Download and prepare the 'genes' data frame.
    attributes <- c("ensembl_gene_id", "ensembl_transcript_id")
    bm_table <- getBM(attributes, filters=filters, values=values, mart=mart)
    genes_tx_id <- as.integer(factor(bm_table$ensembl_transcript_id,
                                     levels=transcripts_tx_name))
    genes <- data.frame(
        tx_id=genes_tx_id,
        gene_id=bm_table$ensembl_gene_id
    )
    ## Call makeTranscriptDb().
    makeTranscriptDb(transcripts, splicings, genes)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method for TranscriptDb objects.
###

setMethod("show", "TranscriptDb",
    function(object)
    {
        cat("TranscriptDb object:\n")
        metadata <- dbReadTable(object@conn, "metadata")
        creation_date <- metadata[metadata$name == "CreationDate", "value"]
        transcript_nrow <- metadata[metadata$name == "transcript_nrow", "value"]
        exon_nrow <- metadata[metadata$name == "exon_nrow", "value"]
        cds_nrow <- metadata[metadata$name == "cds_nrow", "value"]
        cat("| Creation Date: ", creation_date, "\n", sep="")
        cat("| Nb of rows in transcript table: ", transcript_nrow, "\n", sep="")
        cat("| Nb of rows in exon table: ", exon_nrow, "\n", sep="")
        cat("| Nb of rows in cds table: ", cds_nrow, "\n", sep="")
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
        ## Retrieve the "transcripts" element.
        sql <- paste(
            "SELECT transcript._tx_id AS tx_id, tx_name,",
            "tx_chrom, tx_strand, tx_start, tx_end",
            "FROM transcript INNER JOIN transcript_rtree",
            "ON (transcript._tx_id=transcript_rtree._tx_id)",
            "ORDER BY tx_chrom, tx_strand, tx_start, tx_end, tx_id"
        )
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
            "FROM transcript INNER JOIN transcript_rtree",
            "ON (transcript._tx_id=transcript_rtree._tx_id)",
            "INNER JOIN splicing",
            "ON (transcript._tx_id=splicing._tx_id)",
            "INNER JOIN exon",
            "ON (splicing._exon_id=exon._exon_id)",
            "INNER JOIN exon_rtree",
            "ON (exon._exon_id=exon_rtree._exon_id)",
            "LEFT JOIN cds",
            "ON (splicing._cds_id=cds._cds_id)",
            "LEFT JOIN cds_rtree",
            "ON (cds._cds_id=cds_rtree._cds_id)",
            "ORDER BY tx_chrom, tx_strand, tx_start, tx_end, tx_id, exon_rank"
        )
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
            "FROM transcript INNER JOIN transcript_rtree",
            "ON (transcript._tx_id=transcript_rtree._tx_id)",
            "INNER JOIN gene",
            "ON (transcript._tx_id=gene._tx_id)",
            "ORDER BY tx_chrom, tx_strand, tx_start, tx_end, tx_id, gene_id"
        )
        genes <- dbGetQuery(x@conn, sql)
        COL2CLASS <- c(
             tx_id="integer",
             gene_id="character"
        )
        genes <- setDataFrameColClass(genes, COL2CLASS)

        list(transcripts=transcripts, splicings=splicings, genes=genes)
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

