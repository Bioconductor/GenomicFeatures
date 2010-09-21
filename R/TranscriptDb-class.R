### =========================================================================
### TranscriptDb objects
### -------------------------------------------------------------------------


### Not exported.
DB_TYPE_NAME <- "Db type"
DB_TYPE_VALUE <- "TranscriptDb"  # same as the name of the class

### Not exported.
makeFeatureColnames <- function(feature_shortname)
{
    suffixes <- c("_id", "_name", "_chrom", "_strand", "_start", "_end")
    prefixes <- c("_", rep.int("", length(suffixes) - 1L))
    paste(prefixes, feature_shortname, suffixes, sep="")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level accessor (not exported).
###

.getConn <- function(envir) get("conn", envir=envir, inherits=FALSE)

txdbConn <- function(txdb) .getConn(txdb@envir)


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
    data0 <- dbEasyQuery(conn, sql0)
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
                 " WHERE name = '", DB_TYPE_NAME, "'", sep="")
    data <- dbEasyQuery(conn, sql)
    if (nrow(data) != 1L) {
        msg <- paste("the metadata table in the DB has 0 or more ",
                     "than 1 '", DB_TYPE_NAME, "' entries", sep="")
        return(msg)
    }
    db_type <- data[["value"]]
    if (is.na(db_type) || db_type != DB_TYPE_VALUE) {
        msg <- paste("'", DB_TYPE_NAME, "' is not \"", DB_TYPE_VALUE,
                     "\"", sep="")
        return(msg)
    }
    NULL
}

### TODO: Add more checks!
.valid.transcript.table <- function(conn)
{
    colnames <- makeFeatureColnames("tx")
    msg <- .valid.table.colnames(conn, "transcript", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

### TODO: Add more checks!
.valid.exon.table <- function(conn)
{
    colnames <- makeFeatureColnames("exon")
    msg <- .valid.table.colnames(conn, "exon", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

### TODO: Add more checks!
.valid.cds.table <- function(conn)
{
    colnames <- makeFeatureColnames("cds")
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
    colnames <- c("_chrom_id", "chrom", "length", "is_circular")
    msg <- .valid.table.colnames(conn, "chrominfo", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

.valid.TranscriptDb <- function(x)
{
    conn <- txdbConn(x)
    c(.valid.metadata.table(conn),
      .valid.transcript.table(conn),
      .valid.exon.table(conn),
      .valid.cds.table(conn),
      .valid.splicing.table(conn),
      .valid.gene.table(conn),
      .valid.chrominfo.table(conn))
}

setValidity2("TranscriptDb", .valid.TranscriptDb)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor (not exported).
###

TranscriptDb <- function(conn)
{
    if (!is(conn, "SQLiteConnection"))
        stop("'conn' must be an SQLiteConnection object")
    envir <- new.env(parent=emptyenv())
    assign("conn", conn, envir=envir)
    reg.finalizer(envir, function(e) dbDisconnect(.getConn(e)))
    new("TranscriptDb", envir=envir)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Saving/loading.
###

saveFeatures <- function(x, file)
{
    if (!is(x, "TranscriptDb"))
        stop("'x' must be a TranscriptDb object")
    if (!isSingleString(file))
        stop("'file' must be a single string")
    sqliteCopyDatabase(txdbConn(x), file)
}

loadFeatures <- function(file)
{
    if (!isSingleString(file))
        stop("'file' must be a single string")
    if(!file.exists(file))
        stop("file '", file, "' does not exist")
    conn <- dbConnect(SQLite(), file)
    TranscriptDb(conn)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors.
###

setMethod("metadata", "TranscriptDb",
    function(x) dbReadTable(txdbConn(x), "metadata")
)

.getChromInfo <- function(x)
{
    sql <- "SELECT chrom, length, is_circular FROM chrominfo ORDER BY _chrom_id"
    chrominfo <- dbEasyQuery(txdbConn(x), sql)
    COL2CLASS <- c(
         chrom="character",
         length="integer",
         is_circular="logical"
    )
    setDataFrameColClass(chrominfo, COL2CLASS)
}

setMethod("seqnames", "TranscriptDb",
    function(x)
    {
        data <- .getChromInfo(x)
        as.character(data[["chrom"]])
    }
)

setMethod("seqlengths", "TranscriptDb",
    function(x)
    {
        data <- .getChromInfo(x)
        ans <- as.integer(data[["length"]])
        names(ans) <- as.character(data[["chrom"]])
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method for TranscriptDb objects.
###

setMethod("show", "TranscriptDb",
    function(object)
    {
        cat("TranscriptDb object:\n")
        metadata <- metadata(object)
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
        transcripts <- dbEasyQuery(txdbConn(x), sql)
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
        splicings <- dbEasyQuery(txdbConn(x), sql)
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
        genes <- dbEasyQuery(txdbConn(x), sql)
        COL2CLASS <- c(
             tx_id="integer",
             gene_id="character"
        )
        genes <- setDataFrameColClass(genes, COL2CLASS)

        ## Retrieve the "chrominfo" element.
        chrominfo <- .getChromInfo(x)

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

