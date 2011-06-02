### =========================================================================
### TranscriptDb objects
### -------------------------------------------------------------------------


### Not exported.
DB_TYPE_NAME <- "Db type"
DB_TYPE_VALUE <- "TranscriptDb"  # same as the name of the class
DB_SCHEMA_VERSION <- "1.0"
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

.getMetaValue <- function(conn, name) {
    colnames <- c("name", "value")
    msg <- .valid.table.colnames(conn, "metadata", colnames)
    if (!is.null(msg))
        stop(msg)
    sql <-  paste("SELECT * FROM metadata",
                 " WHERE name = '", name, "'", sep="")
    data <- dbEasyQuery(conn, sql)
    if (nrow(data) != 1L) {
        msg <- paste("The metadata table in the DB has 0 or more ",
                     "than 1 '", name, "' entries", sep="")
        stop(msg)
    }
    data$value
}


.valid.colnames <- function(conn, tablename, colnames)
{
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

.valid.table.colnames <- function(conn, tablename, colnames)
{
    tmp <- try(dbExistsTable(conn, tablename), silent=TRUE)
    if (is(tmp, "try-error"))
        return("invalid DB file")
    if (!tmp)
        return(paste("the DB has no ", tablename, " table", sep=""))
    .valid.colnames(conn, tablename, colnames)
}



.valid.metadata.table <- function(conn, type)
{
    colnames <- c("name", "value")
    msg <- .valid.table.colnames(conn, "metadata", colnames)
    if (!is.null(msg))
        return(msg)
    db_type <- try(.getMetaValue(conn, DB_TYPE_NAME), silent = TRUE)
    if(is(db_type, "try-error"))
        return(db_type[1])
    if (is.na(db_type) || db_type != type) {
        msg <- paste("'", DB_TYPE_NAME, "' is not \"", type,
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
    c(.valid.metadata.table(conn, DB_TYPE_VALUE),
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
    seqNames <- .getChromInfo(conn)$chrom
    seqNVals <- rep(TRUE, length(seqNames))
    names(seqNVals) <- seqNames
    new("TranscriptDb", envir=envir, activeSeqs=seqNVals)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Saving/loading. 
###

setMethod("saveFeatures", "TranscriptDb",
          function(x, file)
          {
            if (!is(x, "TranscriptDb"))
              stop("'x' must be a TranscriptDb object")
            if (!isSingleString(file))
              stop("'file' must be a single string")
            sqliteCopyDatabase(txdbConn(x), file)
          }
)

loadFeatures <- function(file)
{
    if (!isSingleString(file))
        stop("'file' must be a single string")
    if(!file.exists(file))
        stop("file '", file, "' does not exist")

    conn <- dbConnect(SQLite(), file)
    if(dbExistsTable(conn, "metadata")) {
        type <- .getMetaValue(conn, "Db type")
            if(type == "TranscriptDb") {
              version <- try(.getMetaValue(conn,"DBSCHEMAVERSION"),
                             silent = TRUE)
              if(is(version, "try-error")){
                  conn <- .fixOldDbSchema(conn)
              }
              return(TranscriptDb(conn))
            }else if(type == "FeatureDb") {
              return(FeatureDb(conn))
            }else{
              stop("The file you are trying to load is of unknown Db type")
            }
    }
    ##TranscriptDb(conn)
}



.fixOldDbSchema <- function(conn) {
    db <- dbConnect(SQLite(), dbname = ":memory:")
    sqliteCopyDatabase(conn, db)
    dbDisconnect(conn)
    sql <- "SELECT  * from chrominfo"
    chromInfo <- dbEasyQuery(db, sql)
    nr <- nrow(chromInfo)
    if( !"is_circular" %in% colnames(chromInfo)){
        is_circular <- rep(NA, nr)
        sql <- paste("ALTER TABLE chrominfo ADD is_circular", is_circular, "INTEGER", sep = " ")
        dbEasyQuery(db, sql)
        sql <- paste("INSERT INTO metadata VALUES('DBSCHEMAVERSION', '",
                DB_SCHEMA_VERSION,"')", sep = "")
        dbEasyQuery(db, sql)
        message("The TranscriptDb object has been updated to the latest schema version 1.0.")
        message("The updated object can be saved using saveFeatures()")
    } 
    db
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
    chrominfo <- dbEasyQuery(x, sql)
    COL2CLASS <- c(
         chrom="character",
         length="integer",
         is_circular="logical"
    )
    setDataFrameColClass(chrominfo, COL2CLASS)
}

setMethod("seqinfo", "TranscriptDb",
    function(x)
    {
        data <- .getChromInfo(txdbConn(x))
        Seqinfo(seqnames = data[["chrom"]],
                seqlengths = data[["length"]],
                isCircular = data[["is_circular"]])
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
        chrominfo <- .getChromInfo(txdbConn(x))

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


### =========================================================================
### Setters and getters for activeSeqs
### -------------------------------------------------------------------------
setGeneric("activeSeqs", function(x) standardGeneric("activeSeqs"))
setGeneric("activeSeqs<-",function(x, value) standardGeneric("activeSeqs<-"))

## setters 
.setSeqNames <- function(x, value){
  ## Must check to make sure that the values are legitimate
  seqNames <- seqlevels(x)
  ## The names must be all containe in seqNames
  if(length(intersect(names(value),seqNames)) == length(value) &&
     ##length(value) == length(seqNames) && ## cannot be shorter than seqNames
     is.logical(value)){ ## and it must be a logical
    x@activeSeqs[names(value)] <- value	
  }else{stop("The replacement value for activeSeqs must be a logical ",
             "vector, with names that match the seqlevels of the ",
             "TranscriptDb object.")
  }
  x
}

setReplaceMethod("activeSeqs","TranscriptDb",
	  function(x, value){.setSeqNames(x,value)})


## getters
setMethod("activeSeqs", "TranscriptDb", function(x){x@activeSeqs})


## TODO: make manual pages.

## library(GenomicFeatures);example(loadFeatures); keep = "chr1"; activeSeqs(txdb)[c("chr1", "chr3")] <- FALSE ; activeSeqs(txdb)

## activeSeqs(txdb)

## library(GenomicFeatures);example(loadFeatures); keep = "chr1"

## How to use the accessors:  (add to manual page)
## activeSeqs(txdb)[value] <- TRUE
## activeSeqs(txdb)[c("chr1", "chr3")] <- FALSE

## activeSeqs(txdb)[seqlevels(txdb)] <- FALSE
## activeSeqs(txdb)[c("chr1", "chr3")] <- TRUE

## then transcripts and transcriptsBy etc will pay attention to non-default vals
