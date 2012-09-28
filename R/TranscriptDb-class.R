### =========================================================================
### TranscriptDb objects
### -------------------------------------------------------------------------


## This is to try and tidy up before setRefClass()
gc()

### Concrete GenomicFeatures types
.TranscriptDb <-
    setRefClass("TranscriptDb", contains="AnnotationDb",
        fields=list(isActiveSeq="logical", seqnameStyle="character"),
        methods=list(
          initialize=function(...) {
              callSuper(...)
              .self$seqnameStyle <- character()
              if (0L == length(dbListTables(conn))) {
                  .self$isActiveSeq <- logical()
              } else {
                  seqNames <- load_chrominfo(.self, set.col.class=TRUE)$chrom
                  .self$isActiveSeq <-
                      structure(!logical(length(seqNames)), .Names=seqNames)
              }
          .self
      }))

### Not exported.
DB_TYPE_NAME <- "Db type"
DB_TYPE_VALUE <- "TranscriptDb"  # same as the name of the class
DB_SCHEMA_VERSION <- "1.0"


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level data loaders.
###
### For internal use only (i.e. NOT exported)
###

load_chrominfo <- function(txdb, set.col.class=FALSE)
{
    sql <- c("SELECT chrom, length, is_circular",
             "FROM chrominfo ORDER BY _chrom_id")
    ans <- queryAnnotationDb(txdb, sql)
    if (!set.col.class)
        return(ans)
    COL2CLASS <- c(
         chrom="character",
         length="integer",
         is_circular="logical"
    )
    setDataFrameColClass(ans, COL2CLASS)
}

load_transcripts <- function(txdb, set.col.class=FALSE)
{
    sql <- c("SELECT _tx_id AS tx_id, tx_name,",
             "  tx_chrom, tx_strand, tx_start, tx_end",
             "FROM transcript",
             "ORDER BY tx_id")
    ans <- queryAnnotationDb(txdb, sql)
    if (!set.col.class)
        return(ans)
    COL2CLASS <- c(
        tx_id="integer",
        tx_name="character",
        tx_chrom="factor",
        tx_strand="factor",
        tx_start="integer",
        tx_end="integer"
    )
    setDataFrameColClass(ans, COL2CLASS)
}

load_splicings <- function(txdb, set.col.class=FALSE)
{
    sql <- c("SELECT _tx_id AS tx_id, exon_rank,",
             "  splicing._exon_id AS exon_id, exon_name,",
             "  exon_chrom, exon_strand, exon_start, exon_end,",
             #"  splicing._cds_id AS cds_id, cds_name,",
             "  splicing._cds_id AS cds_id,",
             "  cds_start, cds_end",
             "FROM splicing",
             "  INNER JOIN exon",
             "    ON (exon_id=exon._exon_id)",
             "  LEFT JOIN cds",
             "    ON (cds_id=cds._cds_id)",
             "ORDER BY tx_id, exon_rank")
    ans <- queryAnnotationDb(txdb, sql)
    if (!set.col.class)
        return(ans)
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
    setDataFrameColClass(ans, COL2CLASS)
}

load_genes <- function(txdb, set.col.class=FALSE)
{
    sql <- c("SELECT transcript._tx_id AS tx_id, gene_id",
             "FROM transcript",
             "  INNER JOIN gene",
             "    ON (transcript._tx_id=gene._tx_id)",
             "ORDER BY tx_id, gene_id")
    ans <- queryAnnotationDb(txdb, sql)
    if (!set.col.class)
        return(ans)
    COL2CLASS <- c(
        tx_id="integer",
        gene_id="character"
    )
    setDataFrameColClass(ans, COL2CLASS)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity of a TranscriptDb object.
###

### Not exported.
makeFeatureColnames <- function(feature_shortname)
{
    suffixes <- c("_id", "_name", "_chrom", "_strand", "_start", "_end")
    prefixes <- c("_", rep.int("", length(suffixes) - 1L))
    paste0(prefixes, feature_shortname, suffixes)
}

### The table specified in a call to .valid.table.colnames() must have at
### least the cols specified in 'colnames'. It's OK if it has more cols or
### if it has them in a different order.

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

    c(AnnotationDbi:::.valid.metadata.table(conn, DB_TYPE_NAME,
                                            DB_TYPE_VALUE),
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
### Accessors.
###

.seqinfo.TranscriptDb <- function(x)
{
    data <- load_chrominfo(x, set.col.class=TRUE)
    ans <- Seqinfo(seqnames=data[["chrom"]],
                   seqlengths=data[["length"]],
                   isCircular=data[["is_circular"]])
    ## also get the genome information.
    sql <- "SELECT value FROM metadata WHERE name='Genome'"
    genome <- unlist(queryAnnotationDb(x, sql))
    names(genome) <- NULL
    if (length(genome) != 0L)
        genome(ans) <- genome
    ans
}

setMethod("seqinfo", "TranscriptDb", .seqinfo.TranscriptDb)

### Setters and getters for isActiveSeq:

### generic internal setter
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

setGeneric("isActiveSeq", function(x) standardGeneric("isActiveSeq"))

setMethod("isActiveSeq", "TranscriptDb", function(x){x$isActiveSeq})

setGeneric("isActiveSeq<-",function(x, value) standardGeneric("isActiveSeq<-"))

setReplaceMethod("isActiveSeq","TranscriptDb",
	  function(x, value){.setSeqNames(x,value)})


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
        transcripts <- load_transcripts(x, set.col.class=TRUE)
        splicings <- load_splicings(x, set.col.class=TRUE)
        genes <- load_genes(x, set.col.class=TRUE)
        chrominfo <- load_chrominfo(x, set.col.class=TRUE)
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
  mcols(bed)$thick <- unlist(ranges(cds_tx), use.names=FALSE)
  bed
})

setMethod("asGFF", "TranscriptDb", function(x) {
  tx_gene <- transcriptsBy(x)
  gene <- unlist(range(tx_gene))
  mcols(gene)$Parent <- CharacterList(character())
  addPrefix <- function(ids, prefix) {
    if (is(ids, "List"))
      ids <- as(ids, "CharacterList")
    prefix <- paste(prefix, ":", sep = "")
    sub(paste("^|^", prefix, sep = ""), prefix, ids)
  }
  mcols(gene)$ID <- addPrefix(names(gene), "GeneID")
  mcols(gene)$Name <- names(gene)
  mcols(gene)$type <- "gene"
  tx <- transcripts(x, columns = c(Parent = "gene_id", ID = "tx_id",
                         Name = "tx_name"))
  mcols(tx)$Parent <- addPrefix(mcols(tx)$Parent, "GeneID")
  mcols(tx)$ID <- addPrefix(mcols(tx)$ID, "TxID")
  mcols(tx)$type <- "mRNA"
  exon <- exons(x, columns = c(Parent = "tx_id"))
  mcols(exon)$Parent <- addPrefix(mcols(exon)$Parent, "TxID")
  mcols(exon)$ID <- NA
  mcols(exon)$Name <- NA
  mcols(exon)$type <- "exon"
  mcols(exon)$Parent <- as(mcols(exon)$Parent, "CharacterList")
  cds <- cds(x, columns = c(Parent = "tx_id"))
  mcols(cds)$Parent <- addPrefix(mcols(cds)$Parent, "TxID")
  mcols(cds)$ID <- NA
  mcols(cds)$Name <- NA
  mcols(cds)$type <- "CDS"
  gff <- c(gene, tx, exon, cds)
  names(gff) <- NULL
  gff
})

setMethod(rtracklayer:::bestFileFormat, "TranscriptDb", function(x) "gff3")

