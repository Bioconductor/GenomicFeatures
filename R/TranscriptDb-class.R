### =========================================================================
### TranscriptDb objects
### -------------------------------------------------------------------------


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
## extractTranscriptsFromGenome because these rely on the names coming
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
## seqlevels(txdb) <- c(chr5 = "5")

### And to reset we do this:
## txdb <- restoreSeqlevels(txdb)


## Bugs :
## This works:
## seqlevels(txdb) <- c(chr4 = "4", chr5 = "5", chr6="6")
## But this doesn't: (needs to go based on match internally)
## seqlevels(txdb) <- c(chr5 = "5", chr6="6", chr4="4")



## Some unit tests needed for the following:

## This should fail: add seqlevels (and it does)
## seqlevels(txdb) <- c(foo = "2")
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
.TranscriptDb <-
    setRefClass("TranscriptDb", contains="AnnotationDb",
        fields=list(.chrom="character",
                    isActiveSeq="logical",
                    seqnameStyle="character",
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
                  .self$seqnameStyle <- character()
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

## seqinfo getter needs to rename and re-sort things based on new2old
## integers every single time it is accessed
.seqinfo.TranscriptDb <- function(x)
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

setMethod("seqinfo", "TranscriptDb", function(x){.seqinfo.TranscriptDb(x)})

### This is a restricted "seqinfo<-" method for TranscriptDb objects
### that only supports replacement of the sequence names or dropping
### and resorting of the sequence names.  Since the getter above is
### resorting and renaming every time, the setter only needs to put
### the new2old field together correctly (in the object).

.seqinfo.TranscriptDbReplace <- function(x, new2old=NULL, force=FALSE, value)
{
    if (!is(value, "Seqinfo"))
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    IN_THIS_CONTEXT <- paste0("when replacing the 'seqinfo' ",
                              "of a TranscriptDb object")

    ## Get the current seqinfo
    x_seqinfo <- seqinfo(x)
        
    ## make sure new2old is set up
    if (is.null(new2old)) {
        if (!identical(value, x_seqinfo))
            stop("'new2old' must be specified ", IN_THIS_CONTEXT)
        return(x)
    }

    if(!(length(new2old) <= length(x_seqinfo))){
        stop("The replacement value must be either a 1 to 1 replacement ",
             "or a subset of the original set ",
                 IN_THIS_CONTEXT)
    }
    
    ## just always set the new2old value up if it's here.
    x$new2old <- new2old
    
    ## store the names where we always have
    x$.chrom <- seqnames(value)
    x
}

setReplaceMethod("seqinfo", "TranscriptDb", function(x, new2old, force, value){
    .seqinfo.TranscriptDbReplace(x, new2old=new2old, force=force, value)})


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
setMethod("seqlevels0", "TranscriptDb", 
    function(x) x$initialize()
)

setGeneric("isActiveSeq", function(x) standardGeneric("isActiveSeq"))

setMethod("isActiveSeq", "TranscriptDb",
    function(x)
    {
        ans <- x$isActiveSeq
        names(ans) <- x$.chrom
        ans
    }
)

setGeneric("isActiveSeq<-",
    function(x, value) standardGeneric("isActiveSeq<-")
)

.mk_isActiveSeqReplacementValue <- function(x, value)
{
    if (!is.logical(value) || any(is.na(value)))
        stop("the supplied 'isActiveSeq' must be a logical vector with no NAs")
    x_isActiveSeq <- isActiveSeq(x)
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

setReplaceMethod("isActiveSeq","TranscriptDb",
    function(x, value)
    {
        value <- .mk_isActiveSeqReplacementValue(x, value)
        x$isActiveSeq <- unname(value)
        x
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
