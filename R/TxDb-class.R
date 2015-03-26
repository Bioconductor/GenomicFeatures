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

## planned fix for seqinfo(): (so it can support the force argument).
## 1) Add a character vector slot to hold seqnames and their "altered"
## names. OR add a seqinfo object slot to hold a temporary seqinfo
## data (which would then be respected by all the methods?)
## 2) new2old argument: Note that while we can support renaming to
## whatever, we cannot allow users to change the length of the
## sequences.  So there are still some restrictions.  However: we CAN
## allow users to limit scope like in isActiveSeq (a slot that will
## become redundant with this one, by simply stashing NA into values
## that the user has excluded.  OR: if we take the seqinfo approach,
## then we would store the temp seqinfo. and allow users to mess with
## it however they liked, and the methods would have to post-filter
## based on that.
## 3) no replace "value", as we can't really allow the user to set the
## seqinfo to just anything.  Because it's a database, some values are
## not allowed?  OR: do we want to just have the results filtered so
## that the results are in line with a seqinfo. slot?  - currently we
## just don't do the setter.

## Herve made a good point about all this.  And that was that if I
## start allowing all these changes to seqinfo, then I have to filter
## the results.  How much of a mess will that be?  Do we include
## ranges that overlap with the new boundary?  If so do we truncate
## those?
## There is currently code that translates from the old seqinfo to the
## new names (currently only naming is allowed).

## No matter which solution I choose, we don't want to allow filtering
## based on circularity.

## AND: I also need to get vals to behave better as it currently is
## ignoring even the meager amount of things that we are trying to
## support via seqlevels.  And this is because it gets translated into
## a query without 1st being translated back to whatever the DB is
## using (so things will break if you use a name for vals after
## renaming your stuff).

## Other considerations: the implementation of a seqinfo slot could
## create additional complications for methods similar to
## extractTranscriptSeqs because these rely on the names coming
## back being the same as the ones in the DB.  In fact, anywhere that
## we have code that relys on this it is at risk for breaking.  So if
## we want to go down the route of making these methods the same
## across the board, then we have to consider the ramifications of
## that decision.  Specifically, some code may not behave as expected
## after users change these things.  For the base level accessors, we
## can hide that by translating or by filtering based on a maleable
## seqinfo slot.  But for more compound operations, existing code may
## be expecting that these values will always be within the parameters
## that have been allowed up till now (IOW no dropping of seqlevels,
## etc.)

## So: 1st decision: what do we want to allow users to change?
## 1) seqlevels (names) - the only thing that is allowed now.
## 2) "drop" seqlevels, (force=TRUE)   - This is what we want to add.
## 3) change values of seqlengths? Nope
## 4) change values of genomes? Nope
## 5) change values of circularity? Nope
## 6) add seqlevels? Nope
## 7) reorder seqlevels? Yes

## To just add #2, option 1, IOW a translation slot (not a full blown
## seqinfo object) should be enough to do the trick.


### Changing the names works
## seqlevels(txdb) <- as.character(1:93)

## But we also want to support subsetting:
## library(TxDb.Hsapiens.UCSC.hg19.knownGene); txdb=TxDb.Hsapiens.UCSC.hg19.knownGene; seqinfo(txdb)
## seqlevels(txdb, force=TRUE) <- c(chr5 = "5")

### And to reset we do this:
## txdb <- restoreSeqlevels(txdb)


## Bugs :
## This works:
## seqlevels(txdb, force=TRUE) <- c(chr4 = "4", chr5 = "5", chr6="6")
## But this doesn't: (needs to go based on match internally)
## seqlevels(txdb, force=TRUE) <- c(chr5 = "5", chr6="6", chr4="4")


## Some unit tests needed for the following:

## This should fail: add seqlevels (and it does)
## seqlevels(txdb, force=TRUE) <- c(foo = "2")
## This throws an error, so that's good


## These next three do not *actually* make any change the object (so
## these are actually "safe", but they SHOULD still throw an error)

## This should fail: change circ
## seqinfo(txdb) <- seqinfo()
## foo = seqinfo(txdb)
## foo@is_circular = rep(TRUE, 93)
## seqinfo(txdb, new2old=1:93) <- foo

## This should fail: change genome
## foo = seqinfo(txdb)
## foo@genome = rep("hg18", 93)
## seqinfo(txdb, new2old=1:93) <- foo

## This should fail: change seqlengths
## foo = seqinfo(txdb)
## foo@seqlengths = rep(1000L, 93)
## seqinfo(txdb, new2old=1:93) <- foo


### Concrete GenomicFeatures types
.TxDb <-
    setRefClass("TxDb", contains="AnnotationDb",
        fields=list(.chrom="character",
                    isActiveSeq="logical",
                    seqlevelsStyle="character",
                    new2old="integer"),
        methods=list(
          initialize=function(...) {
              callSuper(...)
              if (length(dbListTables(conn) != 0L)) {
                  chrominfo <- load_chrominfo(.self, set.col.class=TRUE)
                  ## set up initial new2old slot:
##                   n2o <- 
##                   names(n2o) <- chrominfo$chrom
                  .self$new2old <- as.integer(1:length(chrominfo$chrom))
                  ## deprecate .chrom and isActiveSeq
                  .self$.chrom <- chrominfo$chrom 
                  .self$isActiveSeq <- !logical(length(.self$.chrom)) 
                  .self$seqlevelsStyle <- character()
              }
          .self
      }))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level data loaders.
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
    sql <- c("SELECT chrom, length, is_circular",
             "FROM chrominfo ORDER BY _chrom_id")
    chrominfo <- queryAnnotationDb(txdb, sql)
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
        if (!hasCol(transcripts, "tx_name"))
            COL2CLASS <- COL2CLASS[names(COL2CLASS) != "tx_name"]
        if (!hasCol(transcripts, "tx_type"))
            COL2CLASS <- COL2CLASS[names(COL2CLASS) != "tx_type"]
        if (!identical(names(transcripts), names(COL2CLASS)))
            transcripts <- transcripts[names(COL2CLASS)]
        if (set.col.class)
            transcripts <- setDataFrameColClass(transcripts, COL2CLASS)
    }
    if (drop.tx_name && hasCol(transcripts, "tx_name") &&
                        all(is.na(transcripts$tx_name)))
        transcripts$tx_name <- NULL
    if (drop.tx_type && hasCol(transcripts, "tx_type") &&
                        all(is.na(transcripts$tx_type)))
        transcripts$tx_type <- NULL
    transcripts
}

load_transcripts <- function(txdb, drop.tx_name=FALSE,
                                   drop.tx_type=FALSE,
                                   set.col.class=FALSE)
{
    sql <- "SELECT _tx_id AS tx_id, tx_name, "
    has_tx_type <- .schema_version(dbconn(txdb)) == "1.1"
    if (has_tx_type)
        sql <- c(sql, "tx_type, ")
    sql <- c(sql, "tx_chrom, tx_strand, tx_start, tx_end",
             "FROM transcript",
             "ORDER BY tx_id")
    transcripts <- queryAnnotationDb(txdb, sql)
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
        if (!hasCol(splicings, "cds_id"))
            splicings$cds_id <- NA
        if (!hasCol(splicings, "cds_start"))
            splicings$cds_start <- NA
        if (!hasCol(splicings, "cds_end"))
            splicings$cds_end <- NA
        if (!hasCol(splicings, "exon_name"))
            COL2CLASS <- COL2CLASS[names(COL2CLASS) != "exon_name"]
        if (!hasCol(splicings, "cds_name"))
            COL2CLASS <- COL2CLASS[names(COL2CLASS) != "cds_name"]
        if (!identical(names(splicings), names(COL2CLASS)))
            splicings <- splicings[names(COL2CLASS)]
        if (set.col.class)
            splicings <- setDataFrameColClass(splicings, COL2CLASS)
    }
    if (drop.exon_name && hasCol(splicings, "exon_name") &&
                          all(is.na(splicings$exon_name)))
        splicings$exon_name <- NULL
    if (drop.cds_name && hasCol(splicings, "cds_name") &&
                         all(is.na(splicings$cds_name)))
        splicings$cds_name <- NULL
    splicings
}

load_splicings <- function(txdb, drop.exon_name=FALSE,
                                 drop.cds_name=FALSE,
                                 set.col.class=FALSE)
{
    sql <- c("SELECT _tx_id AS tx_id, exon_rank,",
             "  splicing._exon_id AS exon_id, exon_name,",
             "  exon_chrom, exon_strand, exon_start, exon_end,",
             "  splicing._cds_id AS cds_id, cds_name,",
             "  cds_start, cds_end",
             "FROM splicing",
             "  INNER JOIN exon",
             "    ON (exon_id=exon._exon_id)",
             "  LEFT JOIN cds",
             "    ON (cds_id=cds._cds_id)",
             "ORDER BY tx_id, exon_rank")
    splicings <- queryAnnotationDb(txdb, sql)
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
    sql <- c("SELECT transcript._tx_id AS tx_id, gene_id",
             "FROM transcript",
             "  INNER JOIN gene",
             "    ON (transcript._tx_id=gene._tx_id)",
             "ORDER BY tx_id, gene_id")
    genes <- queryAnnotationDb(txdb, sql)
    .format_genes(genes, set.col.class=set.col.class)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity of a TxDb object.
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
    conn <- AnnotationDbi:::dbconn(x)

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
### Low-level constructor (filled out so that we can be consistent with
### AnnnotationDb
###

## Legacy constructor (still not exported) - used in some other places so
## preserved for now
TxDb <- function(conn)
{
    .TxDb$new(conn=conn)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors.
###

### The "species" method currently defined in AnnotationDbi (1.29.20) for
### AnnotationDb objects is too messed up so here we overwrite it for TxDb
### objects and deprecate it in favor of organism(). Once things are fixed
### in AnnotationDbi, we won't need this anymore and can remove it.
setMethod("species", "TxDb",
    function(object)
    {
         msg <- c("  Calling species() on a ", class(object), " object ",
                  "is *deprecated*.\n  Please use organism() instead.")
        .Deprecated(msg=msg)
        organism(object)
    }
)
setMethod("organism", "TxDb",
    function(object)
    {
        metadata <- metadata(object)
        metadata <- setNames(metadata[ , "value"],
                             tolower(metadata[ , "name"]))
        metadata[["organism"]]
    }
)

## seqinfo getter needs to rename and re-sort things based on new2old
## integers every single time it is accessed
.seqinfo.TxDb <- function(x)
{
    data <- load_chrominfo(x, set.col.class=TRUE)
    ## We take the seqnames from x's private field '.chrom'.
    ans <- Seqinfo(seqnames=x$.chrom, ## stored correctly already
                   seqlengths=data[["length"]][x$new2old],
                   isCircular=data[["is_circular"]][x$new2old])
    ## Then re-arrange the ans value based on the new2old values
    ## also get the genome information.
    sql <- "SELECT value FROM metadata WHERE name='Genome'"
    genome <- unlist(queryAnnotationDb(x, sql))
    names(genome) <- NULL
    if (length(genome) != 0L)
        genome(ans) <- genome
    ans
}

setMethod("seqinfo", "TxDb", function(x){.seqinfo.TxDb(x)})

### This is a restricted "seqinfo<-" method for TxDb objects
### that only supports replacement of the sequence names or dropping
### and resorting of the sequence names.  Since the getter above is
### resorting and renaming every time, the setter only needs to put
### the new2old field together correctly (in the object).

.seqinfo.TxDbReplace <- function(x, new2old=NULL, force=FALSE, value)
{
    if (!is(value, "Seqinfo"))
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    IN_THIS_CONTEXT <- paste0("when replacing the 'seqinfo' ",
                              "of a TxDb object")

    ## Get the current seqinfo
    x_seqinfo <- seqinfo(x)
        
    ## make sure new2old is set up
    if (is.null(new2old)) {
        stop("'new2old' must be specified ", IN_THIS_CONTEXT)
        return(x)
    }
    ## length has to be reasonable to move forward
    if(!is.null(new2old) &&  !(length(new2old) <= length(x_seqinfo))){
        stop("The replacement value must be either a 1 to 1 replacement ",
             "or a subset of the original set ",
                 IN_THIS_CONTEXT)
    }
    
    ## if the new value is smaller, then we need a smaller thing for comparison
    if(!is.null(new2old) && length(new2old) < length(x_seqinfo)){
        equiv_seqinfo <- Seqinfo(seqnames=seqnames(value), 
                   seqlengths=x_seqinfo@seqlengths[new2old],
                   isCircular=x_seqinfo@is_circular[new2old],
                   genome=x_seqinfo@genome[new2old])
    }else{
        equiv_seqinfo <- x_seqinfo
    }

    ## no changes to circ allowed
    if(!identical(value@is_circular, equiv_seqinfo@is_circular)){
        stop("No changes are allowed to circularity ", IN_THIS_CONTEXT)
    }
    
    ## no changes to genome allowed
    if(!identical(value@genome, equiv_seqinfo@genome)){
        stop("No changes are allowed to genome ", IN_THIS_CONTEXT)
    }
    
    ## no changes to lengths allowed
    if(!identical(value@seqlengths, equiv_seqinfo@seqlengths)){
        stop("No changes are allowed to seqlengths ", IN_THIS_CONTEXT)
    }
    
    
    if(force == TRUE && !is.null(new2old)){
        ## just always set the new2old value up if it's here.
        x$new2old <- new2old
        ## and we also need to update the isActiveSeq slot
        x$isActiveSeq <- .isActiveSeq(x)[new2old]
        names(x$isActiveSeq) <- NULL
    }
    if(force != TRUE && !is.null(new2old) &&
       length(new2old) < length(x_seqinfo)){
        stop("You need to use force=TRUE if you want to drop seqlevels.")
    }
    
    ## store the names where we always have
    x$.chrom <- seqnames(value)

    ## And return
    x
}

setReplaceMethod("seqinfo", "TxDb", function(x, new2old, force, value){
    .seqinfo.TxDbReplace(x, new2old=new2old, force=force, value)})


## This is the seqlevels() currently in use from GRanges
## I will have to make my own one of these to fix the problem with the
## new2old being wrong (OR I have to solve it in the replacementmethod
## for seqinfo above...

## ### Default "seqlevels<-" method works on any object 'x' with working
## ### "seqinfo" and "seqinfo<-" methods.
## setReplaceMethod("seqlevels", "ANY",
##     function(x, force=FALSE, value)
##     {
##         ## Make the new Seqinfo object.
##         x_seqinfo <- seqinfo(x)
##         seqlevels(x_seqinfo) <- value
##         ## Map the new sequence levels to the old ones.
##         new2old <- getSeqlevelsReplacementMode(value, seqlevels(x))
##         if (identical(new2old, -2L)) {
##             new2old <- match(value, seqlevels(x))
##         } else if (identical(new2old, -1L)) {
##             new2old <- seq_len(length(value))
##         }
##         ## Do the replacement.
##         seqinfo(x, new2old=new2old, force=force) <- x_seqinfo
##         x
##     }
## )






## Reset seqlevels (seqnames) back to original values.
setMethod("seqlevels0", "TxDb", 
    function(x) x$initialize()
)

setGeneric("isActiveSeq", function(x) standardGeneric("isActiveSeq"))

## , msg="isActiveSeq is deprecated for Bioc 2.13 and above. Please see help(seqinfo) for an alternative approach."

.isActiveSeq <- function(x){
    ans <- x$isActiveSeq
    names(ans) <- x$.chrom
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
        x$isActiveSeq <- unname(value)
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Comparing 2 TxDb objects.
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

