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

## .getConn <- function(envir) get("conn", envir=envir, inherits=FALSE)

## txdbConn <- function(txdb) .getConn(txdb@envir)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity of a TranscriptDb object.
###

### The specified table must have *at least* the cols specified in 'colnames'.
### It's OK if it has more cols or if it has them in a different order.

## .getMetaValue is moved to AnnotationDbi along with a few others
## (commented).  Eventually, the plan is to hide these helper methods inside
## of a reference class so that everyone can use them (but they will still be
## invisible to the users)

## .getMetaValue <- function(conn, name) {
##     colnames <- c("name", "value")
##     msg <- .valid.table.colnames(conn, "metadata", colnames)
##     if (!is.null(msg))
##         stop(msg)
##     sql <-  paste("SELECT * FROM metadata",
##                  " WHERE name = '", name, "'", sep="")
##     data <- dbEasyQuery(conn, sql)
##     if (nrow(data) != 1L) {
##         msg <- paste("The metadata table in the DB has 0 or more ",
##                      "than 1 '", name, "' entries", sep="")
##         stop(msg)
##     }
##     data$value
## }


## .valid.colnames <- function(conn, tablename, colnames)
## {
##     sql0 <- paste("SELECT * FROM ", tablename, " LIMIT 0", sep="")
##     data0 <- dbEasyQuery(conn, sql0)
##     colnames0 <- colnames(data0)
##     if (!all(colnames %in% colnames0)) {
##         msg <- paste("the ", tablename, " table in the DB doesn't have ",
##                      "all the expected columns (",
##                      paste("\"", colnames, "\"", sep="", collapse=", "),
##                      ")", sep="")
##         return(msg)
##     }
##     NULL
## }

## .valid.table.colnames <- function(conn, tablename, colnames)
## {
##     tmp <- try(dbExistsTable(conn, tablename), silent=TRUE)
##     if (is(tmp, "try-error"))
##         return("invalid DB file")
##     if (!tmp)
##         return(paste("the DB has no ", tablename, " table", sep=""))
##     .valid.colnames(conn, tablename, colnames)
## }



## .valid.metadata.table <- function(conn, type)
## {
##     colnames <- c("name", "value")
##     msg <- .valid.table.colnames(conn, "metadata", colnames)
##     if (!is.null(msg))
##         return(msg)
##     db_type <- try(.getMetaValue(conn, DB_TYPE_NAME), silent = TRUE)
##     if(is(db_type, "try-error"))
##         return(db_type[1])
##     if (is.na(db_type) || db_type != type) {
##         msg <- paste("'", DB_TYPE_NAME, "' is not \"", type,
##                      "\"", sep="")
##         return(msg)
##     }
##     NULL
## }

### TODO: Add more checks!
.valid.transcript.table <- function(conn)
{
    colnames <- makeFeatureColnames("tx")
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "transcript", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

### TODO: Add more checks!
.valid.exon.table <- function(conn)
{
    colnames <- makeFeatureColnames("exon")
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "exon", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

### TODO: Add more checks!
.valid.cds.table <- function(conn)
{
    colnames <- makeFeatureColnames("cds")
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "cds", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

### TODO: Add more checks!
.valid.splicing.table <- function(conn)
{
    colnames <- c("_tx_id", "exon_rank", "_exon_id", "_cds_id")
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "splicing", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

### TODO: Add more checks!
.valid.gene.table <- function(conn)
{
    colnames <- c("gene_id", "_tx_id")
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "gene", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

### TODO: Add more checks!
.valid.chrominfo.table <- function(conn)
{
    colnames <- c("_chrom_id", "chrom", "length", "is_circular")
    msg <- AnnotationDbi:::.valid.table.colnames(conn, "chrominfo", colnames)
    if (!is.null(msg))
        return(msg)
    NULL
}

.valid.TranscriptDb <- function(x)
{
    conn <- AnnotationDbi:::dbConn(x)
    c(AnnotationDbi:::.valid.metadata.table(conn, DB_TYPE_VALUE),
      .valid.transcript.table(conn),
      .valid.exon.table(conn),
      .valid.cds.table(conn),
      .valid.splicing.table(conn),
      .valid.gene.table(conn),
      .valid.chrominfo.table(conn))
}

setValidity2("TranscriptDb", .valid.TranscriptDb)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor (filled out so that we can be consistent with
### AnnnotationDb
###

## Legacy constructor (still not exported) - used in some other places so
## preserved for now
TranscriptDb <- function(conn)
{
    .TranscriptDb$new(conn=conn)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Saving/loading. 
###

setMethod("saveFeatures", "TranscriptDb",
          function(x, file)
          {	      
	    .Deprecated(new="saveDb", package="AnnotationDbi")
            saveDb(x, file)
          }
)



## Old loadFeatures will be deprecated, but for now lets just not make any
## current users too unhappy.
loadFeatures <- function(file)
{
    .Deprecated(new="loadDb", package="AnnotationDbi")
    loadDb(file)
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

## setMethod("metadata", "TranscriptDb",
##     function(x) dbReadTable(txdbConn(x), "metadata")
## )


.getChromInfo <- function(x)
{
    sql <- "SELECT chrom, length, is_circular FROM chrominfo ORDER BY _chrom_id"
    chrominfo <- AnnotationDbi:::dbEasyQuery(x, sql)
    COL2CLASS <- c(
         chrom="character",
         length="integer",
         is_circular="logical"
    )
    setDataFrameColClass(chrominfo, COL2CLASS)
}

.seqInfo <- function(x)
{
    data <- .getChromInfo(AnnotationDbi:::dbConn(x))
    ## also get the genome information.
    sql <- "SELECT value FROM metadata WHERE name='Genome'"
    genome <- unlist(AnnotationDbi:::dbEasyQuery(AnnotationDbi:::dbConn(x),sql))
    names(genome) <- NULL
    if(length(genome) > 0){
      Seqinfo(seqnames = data[["chrom"]],
              seqlengths = data[["length"]],
              isCircular = data[["is_circular"]],
              genome = rep(genome,times=length(data[["chrom"]])))
    }else{
      Seqinfo(seqnames = data[["chrom"]],
              seqlengths = data[["length"]],
              isCircular = data[["is_circular"]])
    }
}
setMethod("seqinfo", "TranscriptDb",function(x){.seqInfo(x)})

## library("TxDb.Hsapiens.UCSC.hg19.knownGene");exp = "hg19"; txdb = TxDb.Hsapiens.UCSC.hg19.knownGene; debug(GenomicFeatures:::.seqInfo); seqinfo(txdb)
## checkIdentical(exp, unique(genome(seqinfo(exons(txdb)))))



## Show should just get inherited now.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method for TranscriptDb objects.
###

## setMethod("show", "TranscriptDb",
##     function(object)
##     {
##         cat("TranscriptDb object:\n")
##         metadata <- metadata(object)
##         for (i in seq_len(nrow(metadata))) {
##             cat("| ", metadata[i, "name"], ": ", metadata[i, "value"],
##                 "\n", sep="")
##         }
##     }
## )


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
        transcripts <- AnnotationDbi:::dbEasyQuery(AnnotationDbi:::dbConn(x),
                                                   sql)
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
        splicings <- AnnotationDbi:::dbEasyQuery(AnnotationDbi:::dbConn(x),
                                                 sql)
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
        genes <- AnnotationDbi:::dbEasyQuery(AnnotationDbi:::dbConn(x), sql)
        COL2CLASS <- c(
             tx_id="integer",
             gene_id="character"
        )
        genes <- setDataFrameColClass(genes, COL2CLASS)

        ## Retrieve the "chrominfo" element.
        chrominfo <- .getChromInfo(AnnotationDbi:::dbConn(x))

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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setMethod("asBED", "TranscriptDb", function(x) {
  exons_tx <- exonsBy(x)
  cds_tx <- range(cdsBy(x))
  exons_tx <- exons_tx[names(cds_tx)]
  bed <- asBED(exons_tx)
  values(bed)$thick <- unlist(ranges(cds_tx), use.names=FALSE)
  bed
})

setMethod("asGFF", "TranscriptDb", function(x) {
  tx_gene <- transcriptsBy(x)
  gene <- unlist(range(tx_gene))
  values(gene)$Parent <- CharacterList(character())
  addPrefix <- function(ids, prefix) {
    if (is(ids, "List"))
      ids <- as(ids, "CharacterList")
    prefix <- paste(prefix, ":", sep = "")
    sub(paste("^|^", prefix, sep = ""), prefix, ids)
  }
  values(gene)$ID <- addPrefix(names(gene), "GeneID")
  values(gene)$Name <- names(gene)
  values(gene)$type <- "gene"
  tx <- transcripts(x, columns = c(Parent = "gene_id", ID = "tx_id",
                         Name = "tx_name"))
  values(tx)$Parent <- addPrefix(values(tx)$Parent, "GeneID")
  values(tx)$ID <- addPrefix(values(tx)$ID, "TxID")
  values(tx)$type <- "mRNA"
  exon <- exons(x, columns = c(Parent = "tx_id"))
  values(exon)$Parent <- addPrefix(values(exon)$Parent, "TxID")
  values(exon)$ID <- NA
  values(exon)$Name <- NA
  values(exon)$type <- "exon"
  values(exon)$Parent <- as(values(exon)$Parent, "CharacterList")
  cds <- cds(x, columns = c(Parent = "tx_id"))
  values(cds)$Parent <- addPrefix(values(cds)$Parent, "TxID")
  values(cds)$ID <- NA
  values(cds)$Name <- NA
  values(cds)$type <- "CDS"
  gff <- c(gene, tx, exon, cds)
  names(gff) <- NULL
  gff
})

setMethod(rtracklayer:::bestFileFormat, "TranscriptDb", function(x) "gff3")

### =========================================================================
### Setters and getters for isActiveSeq
### -------------------------------------------------------------------------
## generic internal setter
.setSeqNames <- function(x, value){
  ## Must check to make sure that the values are legitimate
  seqNames <- seqlevels(x)
  ## The names must be all containe in seqNames
  if(length(intersect(names(value),seqNames)) == length(value) &&
     ##length(value) == length(seqNames) && ## cannot be shorter than seqNames
     is.logical(value)){ ## and it must be a logical
    x$isActiveSeq[names(value)] <- value	
  }else{stop("The replacement value for isActiveSeq must be a logical ",
             "vector, with names that match the seqlevels of the object")
  }
  x
}

## setter
setReplaceMethod("isActiveSeq","TranscriptDb",
	  function(x, value){.setSeqNames(x,value)})

## getter
setMethod("isActiveSeq", "TranscriptDb", function(x){x$isActiveSeq})


## fl = "UCSC_knownGene_sample.sqlite"
## library(GenomicFeatures)
## foo = loadDb(fl)
## library(RSQLite)
## conn = dbConnect(SQLite(), fl)
## bar = TranscriptDb(conn)
