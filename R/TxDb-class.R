### =========================================================================
### TxDb objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TxDb schema
###

### Not exported.
DB_TYPE_NAME <- "Db type"
DB_TYPE_VALUE <- "TxDb"  # same as the name of the class below
DB_SCHEMA_VERSION <- "1.1"  # DON'T FORGET TO BUMP THIS WHEN YOU CHANGE THE
                            # SCHEMA

.schema_version <- function(conn)
    numeric_version(AnnotationDbi:::.getMetaValue(conn, "DBSCHEMAVERSION"))

makeFeatureColnames <- function(feature_shortname=c("tx", "exon", "cds"),
                                no.tx_type=FALSE)
{
    feature_shortname <- match.arg(feature_shortname)
    suffixes <- "name"
    if (feature_shortname == "tx" && !no.tx_type)
        suffixes <- c(suffixes, "type")
    suffixes <- c(suffixes, "chrom", "strand", "start", "end")
    ans <- c(paste0("_", feature_shortname, "_id"),
             paste0(feature_shortname, "_", suffixes))
    names(ans) <- c("id", suffixes)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TxDb class definition
###

## This is to try and tidy up before setRefClass()
gc()

### Concrete GenomicFeatures types
.TxDb <-
    setRefClass("TxDb", contains="AnnotationDb",
        fields=list(user_seqlevels="character",
                    user2seqlevels0="integer",
                    isActiveSeq="logical"),
        methods=list(
          initialize=function(...) {
              callSuper(...)
              if (length(dbListTables(conn) != 0L)) {
                  chrominfo <- load_chrominfo(.self, set.col.class=TRUE)
                  .self$user_seqlevels <- chrominfo$chrom 
                  .self$user2seqlevels0 <- seq_along(chrominfo$chrom)
                  ## deprecate isActiveSeq
                  .self$isActiveSeq <- !logical(length(.self$user_seqlevels)) 
              }
          .self
      }))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level data loaders
###
### For internal use only (i.e. NOT exported)
###

.format_chrominfo <- function(chrominfo, set.col.class=FALSE)
{
    COL2CLASS <- c(
        chrom="character",
        length="integer",
        is_circular="logical"
    )
    if (is.null(chrominfo)) {
        chrominfo <- makeZeroRowDataFrame(COL2CLASS)
    } else {
        if (!is.data.frame(chrominfo))
            stop("'chrominfo' must be a data frame")
        if (!identical(names(chrominfo), names(COL2CLASS)))
            chrominfo <- chrominfo[names(COL2CLASS)]
        if (set.col.class)
            chrominfo <- setDataFrameColClass(chrominfo, COL2CLASS)
    }
    chrominfo
}

load_chrominfo <- function(txdb, set.col.class=FALSE)
{
    chrominfo <- TxDb_SELECT_from_chrominfo(txdb)
    colnames(chrominfo) <- sub("^_", "", colnames(chrominfo))
    .format_chrominfo(chrominfo, set.col.class=set.col.class)
}

.format_transcripts <- function(transcripts, drop.tx_name=FALSE,
                                             drop.tx_type=FALSE,
                                             set.col.class=FALSE)
{
    COL2CLASS <- c(
        tx_id="integer",
        tx_name="character",
        tx_type="character",
        tx_chrom="factor",
        tx_strand="factor",
        tx_start="integer",
        tx_end="integer"
    )
    if (is.null(transcripts)) {
        transcripts <- makeZeroRowDataFrame(COL2CLASS)
    } else {
        if (!is.data.frame(transcripts))
            stop("'transcripts' must be a data frame")
        if (!has_col(transcripts, "tx_name"))
            COL2CLASS <- COL2CLASS[names(COL2CLASS) != "tx_name"]
        if (!has_col(transcripts, "tx_type"))
            COL2CLASS <- COL2CLASS[names(COL2CLASS) != "tx_type"]
        if (!identical(names(transcripts), names(COL2CLASS)))
            transcripts <- transcripts[names(COL2CLASS)]
        if (set.col.class)
            transcripts <- setDataFrameColClass(transcripts, COL2CLASS)
    }
    if (drop.tx_name && has_col(transcripts, "tx_name") &&
                        all(is.na(transcripts$tx_name)))
        transcripts$tx_name <- NULL
    if (drop.tx_type && has_col(transcripts, "tx_type") &&
                        all(is.na(transcripts$tx_type)))
        transcripts$tx_type <- NULL
    transcripts
}

load_transcripts <- function(txdb, drop.tx_name=FALSE,
                                   drop.tx_type=FALSE,
                                   set.col.class=FALSE)
{
    transcripts <- TxDb_SELECT_from_transcript(txdb)
    colnames(transcripts) <- sub("^_", "", colnames(transcripts))
    .format_transcripts(transcripts, drop.tx_name=drop.tx_name,
                                     drop.tx_type=drop.tx_type,
                                     set.col.class=set.col.class)
}

.format_splicings <- function(splicings, drop.exon_name=FALSE,
                                         drop.cds_name=FALSE,
                                         set.col.class=FALSE)
{
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
        cds_name="character",
        cds_start="integer",
        cds_end="integer"
    )
    if (is.null(splicings)) {
        splicings <- makeZeroRowDataFrame(COL2CLASS)
    } else {
        if (!is.data.frame(splicings))
            stop("'splicings' must be a data frame")
        if (!has_col(splicings, "cds_id"))
            splicings$cds_id <- NA
        if (!has_col(splicings, "cds_start"))
            splicings$cds_start <- NA
        if (!has_col(splicings, "cds_end"))
            splicings$cds_end <- NA
        if (!has_col(splicings, "exon_name"))
            COL2CLASS <- COL2CLASS[names(COL2CLASS) != "exon_name"]
        if (!has_col(splicings, "cds_name"))
            COL2CLASS <- COL2CLASS[names(COL2CLASS) != "cds_name"]
        if (!identical(names(splicings), names(COL2CLASS)))
            splicings <- splicings[names(COL2CLASS)]
        if (set.col.class)
            splicings <- setDataFrameColClass(splicings, COL2CLASS)
    }
    if (drop.exon_name && has_col(splicings, "exon_name") &&
                          all(is.na(splicings$exon_name)))
        splicings$exon_name <- NULL
    if (drop.cds_name && has_col(splicings, "cds_name") &&
                         all(is.na(splicings$cds_name)))
        splicings$cds_name <- NULL
    splicings
}

load_splicings <- function(txdb, drop.exon_name=FALSE,
                                 drop.cds_name=FALSE,
                                 set.col.class=FALSE)
{
    splicings <- TxDb_SELECT_from_splicings(txdb)
    colnames(splicings) <- sub("^_", "", colnames(splicings))
    .format_splicings(splicings, drop.exon_name=drop.exon_name,
                                 drop.cds_name=drop.cds_name,
                                 set.col.class=set.col.class)
}

.format_genes <- function(genes, set.col.class=FALSE)
{
    COL2CLASS <- c(
        tx_id="integer",
        gene_id="character"
    )
    if (is.null(genes)) {
        genes <- makeZeroRowDataFrame(COL2CLASS)
    } else {
        if (!is.data.frame(genes))
            stop("'genes' must be a data frame")
        if (!identical(names(genes), names(COL2CLASS)))
            genes <- genes[names(COL2CLASS)]
        if (set.col.class)
            genes <- setDataFrameColClass(genes, COL2CLASS)
    }
    genes
}

load_genes <- function(txdb, set.col.class=FALSE)
{
    genes <- TxDb_SELECT_from_gene(txdb)
    colnames(genes) <- sub("^_", "", colnames(genes))
    .format_genes(genes, set.col.class=set.col.class)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity of a TxDb object
###

### The table specified in a call to .valid.table.colnames() must have at
### least the cols specified in 'colnames'. It's OK if it has more cols or
### if it has them in a different order.

### TODO: Add more checks!
.valid.transcript.table <- function(conn)
{
    no_tx_type <- .schema_version(conn) < "1.1"
    colnames <- makeFeatureColnames("tx", no_tx_type)
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

.valid.TxDb <- function(x)
{
    conn <- dbconn(x)

    c(AnnotationDbi:::.valid.metadata.table(conn, DB_TYPE_NAME,
                                            DB_TYPE_VALUE),
      .valid.transcript.table(conn),
      .valid.exon.table(conn),
      .valid.cds.table(conn),
      .valid.splicing.table(conn),
      .valid.gene.table(conn),
      .valid.chrominfo.table(conn))
}

setValidity2("TxDb", .valid.TxDb)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level constructor (not exported)
###

### Only used in makeTxDb().
TxDb <- function(conn) .TxDb$new(conn=conn)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### species() and organism() getters
###

### The "species" method currently defined in AnnotationDbi (1.29.20) for
### AnnotationDb objects is too messed up so here we overwrite it for TxDb
### objects and deprecate it in favor of organism(). Once things are fixed
### in AnnotationDbi, we won't need this anymore and can remove it.
setMethod("species", "TxDb",
    function(object)
    {
         msg <- wmsg("Calling species() on a ", class(object), " object ",
                     "is defunct.\n  Please use organism() instead.")
        .Defunct(msg=msg)
    }
)
setMethod("organism", "TxDb",
    function(object)
    {
        metadata <- metadata(object)
        metadata <- setNames(metadata[ , "value"],
                             tolower(metadata[ , "name"]))
        unname(metadata["organism"])
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqinfo() getter and setter
###

### Get the original seqlevels (i.e. extract them directly from the db).
setMethod("seqlevels0", "TxDb",
    function(x)
    {
        sql <- "SELECT chrom FROM chrominfo ORDER BY _chrom_id"
        ans <- queryAnnotationDb(x, sql)[[1L]]
        attr(ans, "seqlevels0") <- TRUE
        ans
    }
)

### Adapted from default "seqlevels<-" method defined in GenomeInfoDb.
### We only support "renaming" and "strict subsetting" modes.
.set_TxDb_seqlevels <-
    function(x, force=FALSE,
             pruning.mode=c("error", "coarse", "fine", "tidy"),
             value)
{
    x_seqlevels0 <- seqlevels0(x)  # "real" seqlevels (from the db)
    if (identical(value, x_seqlevels0))
        return(x$initialize())
    x_seqlevels <- seqlevels(x)
    ## First we compare the user-supplied seqlevels with 'x_seqlevels0' to
    ## detect the situation where the user intention is to subset the "real"
    ## seqlevels.
    mode <- GenomeInfoDb:::getSeqlevelsReplacementMode(value, x_seqlevels0)
    if (mode == -2L) {
        ## "subsetting of the real seqlevels" mode
        x$user_seqlevels <- value
        x$user2seqlevels0 <- match(value, x_seqlevels0)
        return(x)
    }
    ## Then we compare the user-supplied seqlevels with the current user-
    ## defined seqlevels.
    new2old <- GenomeInfoDb:::getSeqlevelsReplacementMode(value, x_seqlevels)
    if (identical(new2old, -3L)) {
        ## "renaming of user-defined seqlevels" mode
        x$user_seqlevels <- value
        return(x)
    }
    if (identical(new2old, -2L) || identical(new2old, -1L)) {
        ## "subsetting of user-defined seqlevels" mode
        new2old <- match(value, x_seqlevels)
    }
    user2seqlevels0 <- x$user2seqlevels0[new2old]
    na_idx <- which(is.na(user2seqlevels0))
    if (length(na_idx) != 0L) {
        user2seqlevels0[na_idx] <- match(value[na_idx], x_seqlevels0)
        if (anyNA(user2seqlevels0))
            stop(wmsg("adding seqlevels to a TxDb object is not supported"))
        any_dup <- anyDuplicated(user2seqlevels0)
        if (any_dup) {
            idx0 <- user2seqlevels0[any_dup]
            in1string <- paste(value[user2seqlevels0 == idx0], collapse=", ")
            seqlevel0 <- x_seqlevels0[idx0]
            stop(wmsg("more than one user-supplied seqlevels ",
                      "(", in1string, ") refer to the same seqlevel ",
                      "stored in the db (", seqlevel0, ")"))
        }
    }
    x$user_seqlevels <- unname(value)
    x$user2seqlevels0 <- user2seqlevels0
    x
}

setReplaceMethod("seqlevels", "TxDb", .set_TxDb_seqlevels)

get_TxDb_seqinfo0 <- function(x)
{
    data <- load_chrominfo(x, set.col.class=TRUE)
    ans <- Seqinfo(seqnames=data[ , "chrom"],
                   seqlengths=data[ , "length"],
                   isCircular=data[ , "is_circular"])
    sql <- "SELECT value FROM metadata WHERE name='Genome'"
    genome <- unlist(queryAnnotationDb(x, sql))
    names(genome) <- NULL
    if (length(genome) != 0L)
        genome(ans) <- genome
    ans
}

.get_TxDb_seqinfo <- function(x)
{
    seqinfo0 <- get_TxDb_seqinfo0(x)
    ans <- seqinfo0[seqlevels(seqinfo0)[x$user2seqlevels0]]
    seqnames(ans) <- x$user_seqlevels
    ans
}

setMethod("seqinfo", "TxDb", .get_TxDb_seqinfo)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isActiveSeq() getter and setter
###

setGeneric("isActiveSeq", function(x) standardGeneric("isActiveSeq"))

## , msg="isActiveSeq is deprecated for Bioc 2.13 and above. Please see help(seqinfo) for an alternative approach."

.isActiveSeq <- function(x){
    ans <- x$isActiveSeq[x$user2seqlevels0]
    names(ans) <- x$user_seqlevels
    ans
}
setMethod("isActiveSeq", "TxDb",
          function(x){
              #.Deprecated("seqlevels", package="GenomicFeatures")
              .isActiveSeq(x)
          })

setGeneric("isActiveSeq<-",
    function(x, value) standardGeneric("isActiveSeq<-")
)

.mk_isActiveSeqReplacementValue <- function(x, value)
{
    if (!is.logical(value) || any(is.na(value)))
        stop("the supplied 'isActiveSeq' must be a logical vector with no NAs")
    x_isActiveSeq <- .isActiveSeq(x)
    current_names <- names(x_isActiveSeq)
    supplied_names <- names(value)
    if (is.null(supplied_names)) {
        if (length(value) != length(x_isActiveSeq))
            stop("when unnamed, the supplied 'isActiveSeq' must ",
                 "have the same length as the current 'isActiveSeq'")
        names(value) <- current_names
        return(value)
    }
    if (any(duplicated(supplied_names)))
        stop("the supplied 'isActiveSeq' has duplicated names")
    idx <- match(supplied_names, current_names)
    if (any(is.na(idx)))
        stop("the names of the supplied 'isActiveSeq' must ",
             "match the names of the current 'isActiveSeq'")
    x_isActiveSeq[idx] <- value
    x_isActiveSeq
}

setReplaceMethod("isActiveSeq","TxDb",
    function(x, value)
    {
        #.Deprecated("seqlevels", package="GenomicFeatures")
        value <- .mk_isActiveSeqReplacementValue(x, value)
        x$isActiveSeq[x$user2seqlevels0] <- unname(value)
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### keep_user_seqlevels_from_TxDb()
###
### Used by the TxDb extractors to correct the seqlevels of GRanges or
### GRangesList object 'x' before it's returned to the user. The correction
### is made according to the "user defined" and "active" seqlevels that are
### currently set on TxDb object 'txdb'.
### 'x' is assumed to have the original TxDb seqlevels on it.

keep_user_seqlevels_from_TxDb <- function(x, txdb)
{
    txdb_seqlevels0 <- seqlevels0(txdb)
    stopifnot(setequal(seqlevels(x), txdb_seqlevels0))
    from_seqlevels <- txdb_seqlevels0[txdb$user2seqlevels0]
    to_seqlevels <- txdb$user_seqlevels
    new_seqlevels <- setNames(to_seqlevels, from_seqlevels)
    new_seqlevels <- new_seqlevels[txdb$isActiveSeq[txdb$user2seqlevels0]]
    seqlevels(x, pruning.mode="coarse") <- new_seqlevels
    x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Comparing 2 TxDb objects
###

### Used in unit tests for makeTxDbFromGRanges().
.format_txdb_dump <- function(transcripts=NULL,
                              splicings=NULL,
                              genes=NULL,
                              chrominfo=NULL)
{
    transcripts <- .format_transcripts(transcripts, drop.tx_name=TRUE,
                                                    drop.tx_type=TRUE,
                                                    set.col.class=TRUE)
    splicings <- .format_splicings(splicings, drop.exon_name=TRUE,
                                              drop.cds_name=TRUE,
                                              set.col.class=TRUE)
    genes <- .format_genes(genes, set.col.class=TRUE)
    chrominfo <- .format_chrominfo(chrominfo, set.col.class=TRUE)
    list(transcripts=transcripts, splicings=splicings,
         genes=genes, chrominfo=chrominfo)
}

### Dump the entire db into a list of data frames, say 'txdb_dump', that can
### then be used to recreate the original db with 'do.call(makeTxDb, txdb_dump)'
### with no loss of information (except possibly for some of the metadata).
### Note that the transcripts are dumped in the same order in all the
### data frames.
setMethod("as.list", "TxDb",
    function(x, ...)
    {
        transcripts <- load_transcripts(x)
        splicings <- load_splicings(x)
        genes <- load_genes(x)
        chrominfo <- load_chrominfo(x)
        .format_txdb_dump(transcripts, splicings, genes, chrominfo)
    }
)

compareTxDbs <- function(txdb1, txdb2)
{
    if (!is(txdb1, "TxDb") || !is(txdb2, "TxDb"))
        stop("'txdb1' and 'txdb2' must be TxDb objects")
    txdb1_dump <- as.list(txdb1)
    txdb2_dump <- as.list(txdb2)
    identical(txdb1_dump, txdb2_dump)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setMethod("asBED", "TxDb", function(x) {
  exons_tx <- exonsBy(x)
  cds_tx <- range(cdsBy(x))
  exons_tx <- exons_tx[names(cds_tx)]
  bed <- asBED(exons_tx)
  mcols(bed)$thick <- unlist(ranges(cds_tx), use.names=FALSE)
  bed
})

setMethod("asGFF", "TxDb", function(x) {
  addPrefix <- function(ids, prefix) {
    if (is(ids, "List"))
      ids <- as(ids, "CharacterList")
    prefix <- paste(prefix, ":", sep = "")
    sub(paste("^|^", prefix, sep = ""), prefix, ids)
  }

  ## Make 'gene' GRanges.
  tx_gene <- transcriptsBy(x)
  gene <- unlist(range(tx_gene))
  mcols(gene)$Parent <- CharacterList(character())
  mcols(gene)$ID <- addPrefix(names(gene), "GeneID")
  mcols(gene)$Name <- names(gene)
  mcols(gene)$type <- "gene"

  ## Make 'tx' GRanges.
  tx <- transcripts(x, columns = c(Parent = "gene_id", ID = "tx_id",
                         Name = "tx_name"))
  mcols(tx)$Parent <- addPrefix(mcols(tx)$Parent, "GeneID")
  mcols(tx)$ID <- addPrefix(mcols(tx)$ID, "TxID")
  mcols(tx)$type <- "mRNA"

  ## Make 'exon' GRanges.
  exon <- exons(x, columns = c(Parent = "tx_id", Name = "exon_name"))
  mcols(exon)$Parent <- addPrefix(mcols(exon)$Parent, "TxID")
  mcols(exon)$ID <- NA
  mcols(exon)$type <- "exon"
  mcols(exon)$Parent <- as(mcols(exon)$Parent, "CharacterList")
  mcols(exon) <- mcols(exon)[ , colnames(mcols(gene))]

  ## Make 'cds' GRanges.
  cds <- cds(x, columns = c(Parent = "tx_id", Name = "cds_name"))
  mcols(cds)$Parent <- addPrefix(mcols(cds)$Parent, "TxID")
  mcols(cds)$ID <- NA
  mcols(cds)$type <- "CDS"
  mcols(cds) <- mcols(cds)[ , colnames(mcols(gene))]

  gff <- c(gene, tx, exon, cds)
  names(gff) <- NULL
  gff
})

setMethod(rtracklayer:::bestFileFormat, "TxDb", function(x) "gff3")

setMethod("show", "TxDb",
    function(object)
    {
        cat(class(object), "object:\n")
        metadata <- metadata(object)
        for (i in seq_len(nrow(metadata))) {
            cat("# ", metadata[i, "name"], ": ", metadata[i, "value"],
                "\n", sep="")
        }
    }
)

