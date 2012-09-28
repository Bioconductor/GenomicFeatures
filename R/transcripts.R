### =========================================================================
### The transcripts(), exons(), cds() and promoters() extractors
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### A high-level representation of the db relational model used for
### generating SQL queries.
###

### Note that we omit the *_start and *_end cols.
.ALLCOLS <- c("gene_id",
              "tx_id", "tx_name", "tx_chrom", "tx_strand",
              "exon_id", "exon_name", "exon_chrom", "exon_strand",
              "cds_id", "cds_name", "cds_chrom", "cds_strand",
              "exon_rank")

.CORETAGS <- c("id", "chrom", "strand", "start", "end")

### THE TRANSCRIPT CENTRIC POINT OF VIEW:
### If we look at the db from a transcript centric point of view, then the
### graph of relations between the tables becomes a tree where the 'transcript'
### table is the root:
###                            transcript
###                             /      \
###                          gene    splicing
###                                  /      \
###                                exon     cds
###
.TRANSCRIPT_CENTRIC_DBDESC <- list(
    CORECOLS=structure(paste("tx", .CORETAGS, sep="_"), names=.CORETAGS),
    ### Each col defined in the db is assigned to the closest table (starting
    ### from the 'transcript' table) where it can be found. The 'transcript'
    ### element must be the 1st in the list:
    COLMAP=list(
        transcript=c("tx_id", makeFeatureColnames("tx")),
        gene="gene_id",
        splicing=c("exon_rank", "exon_id", "_exon_id", "cds_id", "_cds_id"),
        exon=c("exon_name", "exon_chrom", "exon_strand",
               "exon_start", "exon_end"),
        cds=c("cds_name", "cds_chrom", "cds_strand", "cds_start", "cds_end")
    ),
    ### For each table that is not the root or a leaf table in the above tree,
    ### we list the tables that are below it:
    CHILDTABLES=list(
        splicing=c("exon", "cds")
    ),
    ### For each table that is not the root table in the above tree, we specify
    ### the join condition with the parent table:
    JOINS=c(
        gene="transcript._tx_id=gene._tx_id",
        splicing="transcript._tx_id=splicing._tx_id",
        exon="splicing._exon_id=exon._exon_id",
        cds="splicing._cds_id=cds._cds_id"
    )
)

### THE EXON CENTRIC POINT OF VIEW:
### If we look at the db from an exon centric point of view, then the
### graph of relations between the tables becomes a tree where the 'exon'
### table is the root:
###                               exon
###                                |
###                             splicing
###                            /   |    \
###                   transcript  gene  cds
###
.EXON_CENTRIC_DBDESC <- list(
    CORECOLS=structure(paste("exon", .CORETAGS, sep="_"), names=.CORETAGS),
    ### Each col defined in the db is assigned to the closest table (starting
    ### from the 'exon' table) where it can be found. The 'exon' element must
    ### be the 1st in the list:
    COLMAP=list(
        exon=c("exon_id", makeFeatureColnames("exon")),
        splicing=c("tx_id", "_tx_id", "exon_rank", "cds_id", "_cds_id"),
        transcript=c("tx_name", "tx_chrom", "tx_strand", "tx_start", "tx_end"),
        gene="gene_id",
        cds=c("cds_name", "cds_chrom", "cds_strand", "cds_start", "cds_end")
    ),
    ### For each table that is not the root or a leaf table in the above tree,
    ### we list the tables that are below it:
    CHILDTABLES=list(
        splicing=c("transcript", "gene", "cds")
    ),
    ### For each table that is not the root table in the above tree, we specify
    ### the join condition with the parent table:
    JOINS=c(
        splicing="exon._exon_id=splicing._exon_id",
        transcript="splicing._tx_id=transcript._tx_id",
        gene="splicing._tx_id=gene._tx_id",
        cds="splicing._cds_id=cds._cds_id"
    )
)

### THE CDS CENTRIC POINT OF VIEW:
### If we look at the db from a cds centric point of view, then the
### graph of relations between the tables becomes a tree where the 'cds'
### table is the root:
###                               cds
###                                |
###                             splicing
###                            /   |    \
###                   transcript  gene  exon
###
.CDS_CENTRIC_DBDESC <- list(
    CORECOLS=structure(paste("cds", .CORETAGS, sep="_"), names=.CORETAGS),
    ### Each col defined in the db is assigned to the closest table (starting
    ### from the 'cds' table) where it can be found. The 'cds`' element must
    ### be the 1st in the list:
    COLMAP=list(
        cds=c("cds_id", makeFeatureColnames("cds")),
        splicing=c("tx_id", "_tx_id", "exon_rank", "exon_id", "_exon_id"),
        transcript=c("tx_name", "tx_chrom", "tx_strand", "tx_start", "tx_end"),
        gene="gene_id",
        exon=c("exon_name", "exon_chrom", "exon_strand",
               "exon_start", "exon_end")
    ),
    ### For each table that is not the root or a leaf table in the above tree,
    ### we list the tables that are below it:
    CHILDTABLES=list(
        splicing=c("transcript", "gene", "exon")
    ),
    ### For each table that is not the root table in the above tree, we specify
    ### the join condition with the parent table:
    JOINS=c(
        splicing="cds._cds_id=splicing._cds_id",
        transcript="splicing._tx_id=transcript._tx_id",
        gene="splicing._tx_id=gene._tx_id",
        exon="splicing._exon_id=exon._exon_id"
    )
)

.DBDESC <- list(
    transcript=.TRANSCRIPT_CENTRIC_DBDESC,
    exon=.EXON_CENTRIC_DBDESC,
    cds=.CDS_CENTRIC_DBDESC
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some utility functions.
###

.checkargColumns <- function(columns, valid_columns)
{
    if (is.null(columns))
        return()
    if (is.character(columns) && all(columns %in% valid_columns))
        return()
    valid_columns <- paste0("'", valid_columns, "'", collapse = ", ")
    stop("'columns' must be NULL or a character vector ",
         "with values in ", valid_columns)
}

.getClosestTable <- function(root_table, colnames)
{
    COLMAP <- .DBDESC[[root_table]]$COLMAP
    ans <- character(length(colnames))
    ans[] <- NA_character_
    for (tablename in names(COLMAP))
        ans[colnames %in% COLMAP[[tablename]]] <- tablename
    ans
}

.assignColToClosestTable <- function(root_table, colnames)
{
    COLMAP <- .DBDESC[[root_table]]$COLMAP
    lapply(COLMAP, function(cols) intersect(cols, colnames))
}

.asQualifiedColnames <- function(root_table, colnames)
{
    colnames[colnames == "tx_id"] <- "_tx_id"
    colnames[colnames == "exon_id"] <- "_exon_id"
    colnames[colnames == "cds_id"] <- "_cds_id"
    paste(.getClosestTable(root_table, colnames), colnames, sep=".")
}

.joinRootToChildTables <- function(root_table, child_tables)
{
    COLMAP <- .DBDESC[[root_table]]$COLMAP
    CHILDTABLES <- .DBDESC[[root_table]]$CHILDTABLES
    JOINS <- .DBDESC[[root_table]]$JOINS
    all_tables <- names(COLMAP)
    ans <- all_tables[1L]
    for (i in seq_len(length(all_tables))[-1L]) {
        right_table <- all_tables[i]
        right_children <- c(right_table, CHILDTABLES[[right_table]])
        if (length(intersect(child_tables, right_children)) != 0L)
            ans <- paste(ans, "LEFT JOIN", right_table,
                              "ON", JOINS[[right_table]])
    }
    ans
}

### In the case of TranscriptDb objects, the distance between 'root_table'
### and 'child_table' is always <= 2.
### TODO: Revisit this. Not sure it would be guaranteed to work correctly if
### the distance between 'root_table' and 'child_table' was >= 3.
.joinPrimaryKeyToChildTable <- function(root_table, child_table)
{
    COLMAP <- .DBDESC[[root_table]]$COLMAP
    CHILDTABLES <- .DBDESC[[root_table]]$CHILDTABLES
    JOINS <- .DBDESC[[root_table]]$JOINS
    all_tables <- names(COLMAP)
    ans <- ""
    for (i in seq_len(length(all_tables))[-1L]) {
        right_table <- all_tables[i]
        right_children <- c(right_table, CHILDTABLES[[right_table]])
        if (length(intersect(child_table, right_children)) != 0L) {
            if (ans == "") {
                ans <- right_table
                next
            }
            ans <- paste(ans, "LEFT JOIN", right_table,
                              "ON", JOINS[[right_table]])
        }
    }
    ans
}

### convert a named list into an SQL where condition
.sqlWhereIn <- function(vals)
{
    if (length(vals) == 0L)
        return("")
    sql <-
      lapply(seq_len(length(vals)), function(i) {
               v <- vals[[i]]
               if (!is.numeric(v))
                 v <- paste0("'", v, "'")
               v <- paste0("(", paste(v, collapse=","), ")")
               v <- paste0(names(vals)[i], " IN ", v)
               paste0("(", v, ")")
            })
    paste("WHERE", paste(unlist(sql), collapse = " AND "))
}

.extractData <- function(root_table, txdb, what_cols, child_tables, vals,
                         orderby_cols=NULL)
{
    SQL_what <- paste(.asQualifiedColnames(root_table, what_cols),
                      collapse=", ")
    SQL_from <- .joinRootToChildTables(root_table, child_tables)
    SQL_where <- .sqlWhereIn(vals)
    if (length(orderby_cols) == 0L)
        SQL_orderby <- ""
    else
        SQL_orderby <- paste("ORDER BY", paste(orderby_cols, collapse=", "))
    SQL <- paste("SELECT DISTINCT", SQL_what, "FROM", SQL_from,
                 SQL_where, SQL_orderby)
    ans <- queryAnnotationDb(txdb, SQL)
    names(ans) <- what_cols
    ans
}

.getWhereCols <- function(vals)
{
    if (is.null(vals))
        return(character(0))
    ans <- NULL
    if (is.list(vals)) {
        if (length(vals) == 0L)
            return(character(0))
        ans <- names(vals)
    }
    if (is.null(ans))
        stop("'vals' must be NULL or a named list")
    #valid_columns <- setdiff(.ALLCOLS, "exon_rank")
    valid_columns <- .ALLCOLS
    if (!all(ans %in% valid_columns)) {
        valid_columns <- paste0("'", valid_columns, "'", collapse = ", ")
        stop("'vals' must be NULL or a list with names ",
             "in ", valid_columns)
    }
    ans
}

.extractRootData <- function(root_table, txdb, vals, root_columns)
{
    CORECOLS <- .DBDESC[[root_table]]$CORECOLS
    orderby_cols <- CORECOLS[c("chrom", "strand", "start", "end")]
    what_cols <- unique(c(CORECOLS, root_columns))
    where_cols <- .getWhereCols(vals)
    if (is.list(vals))
        names(vals) <- .asQualifiedColnames(root_table, where_cols)
    where_tables <- unique(.getClosestTable(root_table, where_cols))
    .extractData(root_table, txdb, what_cols, where_tables, vals, orderby_cols)
}

.extractDataFromChildTable <- function(root_table, txdb,
                                       primary_key, ids,
                                       child_table, child_columns)
{
    ans_names <- c(primary_key, child_columns)
    primary_key <- paste0("_", primary_key)
    what_cols <- c(primary_key,
                   .asQualifiedColnames(root_table, child_columns))
    SQL_what <- paste(what_cols, collapse=", ")
    SQL_from <- .joinPrimaryKeyToChildTable(root_table, child_table)
    vals <- list(ids)
    names(vals) <- primary_key
    SQL_where <- .sqlWhereIn(vals)
    SQL_orderby <- SQL_what
    SQL <- paste("SELECT DISTINCT", SQL_what, "FROM", SQL_from,
                 SQL_where, "ORDER BY", SQL_orderby)
    ans <- dbEasyQuery(AnnotationDbi:::dbConn(txdb), SQL)
    names(ans) <- ans_names
    ans
}

.extractChildData <- function(root_table, txdb, ids, assigned_columns)
{
    primary_key <- .DBDESC[[root_table]]$CORECOLS["id"]
    ans <- NULL
    all_tables <- names(assigned_columns)
    for (i in seq_len(length(all_tables))[-1L]) {
        child_columns <- assigned_columns[[i]]
        if (length(child_columns) == 0L)
            next
        child_table <- all_tables[i]
        data0 <- .extractDataFromChildTable(root_table, txdb,
                                            primary_key, ids,
                                            child_table, child_columns)
        data <- lapply(data0[ , -1L, drop=FALSE],
                       function(col0)
                       {
                           col <- split(col0, data0[[1L]])
                           col <- col[as.character(ids)]
                           class0 <- class(col0)
                           class <- paste0(toupper(substr(class0, 1L, 1L)),
                                           substr(class0, 2L, nchar(class0)),
                                           "List")
                           get(class)(unname(col))
                       })
        data <- DataFrame(data)
        if (is.null(ans))
            ans <- data
        else
            ans <- c(ans, data)
    }
    ans
}

## make a named list from the metadata data.frame
.makeMetadataList <- function(meta){
    lst <- as.list(meta[,2])
    names(lst) <- meta[,1]
    lst
}

## assign this to the metadata list in relevant object
.assignMetadataList <- function(obj, txdb){
    metadata(obj)[[1]] <- .makeMetadataList(metadata(txdb))
    names(metadata(obj))[[1]] <- "genomeInfo"
    obj
}

.extractFeatureRowsAsGRanges <- function(root_table, txdb, vals, columns)
{
    CORECOLS <- .DBDESC[[root_table]]$CORECOLS
    forbidden_columns <- c(CORECOLS["chrom"], CORECOLS["strand"])
    .checkargColumns(columns, setdiff(.ALLCOLS, forbidden_columns))
    assigned_columns <- .assignColToClosestTable(root_table, columns)
    root_columns <- assigned_columns[[root_table]]
    ## Extract the data from the db.
    root_data <- .extractRootData(root_table, txdb, vals, root_columns)
    child_data <- .extractChildData(root_table, txdb,
                          root_data[[CORECOLS["id"]]], assigned_columns)
    ## Construct the GRanges object and return it.
    ans_seqinfo <- seqinfo(txdb)
    ans_seqnames <- factor(root_data[[CORECOLS["chrom"]]],
                           levels=seqlevels(ans_seqinfo))
    ans_ranges <- IRanges(start=root_data[[CORECOLS["start"]]],
                          end=root_data[[CORECOLS["end"]]])
    ans_strand <- strand(root_data[[CORECOLS["strand"]]])

    activeNames <- names(isActiveSeq(txdb))[isActiveSeq(txdb)]
    seqinfo <- seqinfo(txdb)[activeNames]
    ans <- GRanges(seqnames = ans_seqnames,
                   ranges = ans_ranges,
                   strand = ans_strand,
                   seqinfo = seqinfo)

    ans_values <- c(DataFrame(root_data[root_columns]), child_data)
    if (is.null(names(columns)))
      names(columns) <- columns
    mcols(ans)[names(columns)] <- ans_values[columns]
    .assignMetadataList(ans, txdb)
}


## This is used to create a list from the isActiveSeq slot 
.makeActiveSeqsList <- function(type, txdb){
    actSqs <- isActiveSeq(txdb)
    keepSeqs <- names(actSqs)[actSqs]
    res <- list(keepSeqs)
    names(res) <- type
    res
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The top-level extractors.
###
### Proposal:
###   - rename the 'vals' arg -> 'filter'
###   - rename the 'columns' arg -> 'colnames'
###

setGeneric("transcripts", function(x, ...) standardGeneric("transcripts"))

setMethod("transcripts", "data.frame",
    function(x, vals=NULL, columns=c("tx_id", "tx_name"))
        stop("Please use 'transcripts_deprecated' for older ",
             "data.frame-based transcript metadata.")
)

setMethod("transcripts", "TranscriptDb",
    function(x, vals=NULL, columns=c("tx_id", "tx_name")){
        vals = c(vals, .makeActiveSeqsList("tx_chrom", x))
        .extractFeatureRowsAsGRanges("transcript", x, vals, columns)
      }
)

setGeneric("exons", function(x, ...) standardGeneric("exons"))

setMethod("exons", "data.frame",
    function(x, vals=NULL, columns="exon_id")
        stop("Please use 'exons_deprecated' for older ",
             "data.frame-based transcript metadata.")
)

setMethod("exons", "TranscriptDb",
    function(x, vals=NULL, columns="exon_id"){
        vals = c(vals, .makeActiveSeqsList("exon_chrom", x))
        .extractFeatureRowsAsGRanges("exon", x, vals, columns)
        }
)

setGeneric("cds", function(x, ...) standardGeneric("cds"))

setMethod("cds", "TranscriptDb",
    function(x, vals=NULL, columns="cds_id"){
        vals = c(vals, .makeActiveSeqsList("cds_chrom", x))
        .extractFeatureRowsAsGRanges("cds", x, vals, columns)
        }
)

setGeneric("promoters", signature="x",
    function(x, upstream=2000, downstream=200, ...)
        standardGeneric("promoters"))

setMethod("promoters", "TranscriptDb",
    function(x, upstream=2000, downstream=200, ...)
    {
        gr <- transcripts(x, ...)
        promoters(gr, upstream=upstream, downstream=downstream, ...)
    }
) 

setMethod("promoters", "GenomicRanges",
    function(x, upstream=2000, downstream=200, ...)
    {
        if (!isSingleNumber(upstream))
            stop("'upstream' must be a single integer")
        if (!is.integer(upstream))
            upstream <- as.numeric(upstream)
        if (!isSingleNumber(downstream))
            stop("'downstream' must be a single integer")
        if (!is.integer(downstream))
            downstream <- as.numeric(downstream)
        if (upstream < 0 | downstream < 0)
            stop("'upstream' and 'downstream' must be integers >= 0")
        if (any(strand(x) == "*"))
            warning("'*' ranges were treated as '+'")
        on_plus <- which(strand(x) == "+" | strand(x) == "*")
        on_plus_TSS <- start(x)[on_plus]
        start(x)[on_plus] <- on_plus_TSS - upstream
        end(x)[on_plus] <- on_plus_TSS + downstream - 1L
        on_minus <- which(strand(x) == "-")
        on_minus_TSS <- end(x)[on_minus]
        end(x)[on_minus] <- on_minus_TSS + upstream
        start(x)[on_minus] <- on_minus_TSS - downstream + 1L
        x
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extractors for features in other databases.
###
### This is for extractors that do NOT point to the TranscriptDb proper.
### Such extractors can point to other databases (mirbase.db) OR they can
### point to other FeatureDbs within the same package.

## helpers for microRNAs

## helpers for for translating chromosomes. Now we have to assume a universal
## translator.  It seems that the chroms are in biomaRt style for mirbase.  So
## for biomaRt, return them as is, but for UCSC, add "chr" prefix.
.translateChromsForUCSC <- function(csomes){
  paste0("chr", csomes)
}

.translateChromsForBiomaRt <- function(csomes){
  csomes
}

.syncSeqlevel <- function(txdb, ans){
  isActSeq <- isActiveSeq(txdb)
  n2oNames <- levels(seqnames(ans))
  n2o <- match(seqnames(seqinfo(txdb)), n2oNames)
  seqinfo(ans, new2old=n2o) <- seqinfo(txdb)
  seqlevels(ans, force=TRUE) <- names(isActSeq)[isActSeq]
  ans
}

## main function
.microRNAs <- function(txdb){
  ## get the data about whether or not we have any info.
  con <- AnnotationDbi:::dbConn(txdb)
  bld <- dbGetQuery(con,
           "SELECT value FROM metadata WHERE name='miRBase build ID'")
  src <- DBI:::dbGetQuery(con,
           "SELECT value FROM metadata WHERE name='Data source'")[[1]]

  ## And if not - bail out with message
  if(is.na(bld) || dim(bld)[1]==0){
    stop("this TranscriptDb does not have a miRBase build ID specified")}
  ## now connect to mirbase
  require(mirbase.db) ## strictly required

  ## What I need is the join of mirna with mirna_chromosome_build (via _id),
  ## that is then filtered to only have rows that match the species which goes
  ## with the build.
  
  ## connection
  mcon <- mirbase_dbconn()
  ## 1st lets get the organism abbreviation
  sql <- paste0("SELECT organism FROM mirna_species WHERE genome_assembly='",
                bld, "'")
  organism <- dbGetQuery(mcon, sql)[[1]]
  ## now get data and make a GRanges from it
  sql <- paste0("SELECT * from mirna_chromosome_build AS csome INNER JOIN ",
                "(SELECT _id,mirna_id,organism from mirna) AS mirna ",
                "WHERE mirna._id=csome._id AND organism='", organism, "' ")
  data <- dbGetQuery(mcon, sql)

  ## convert chromosomes
  csomes <- switch(src,
                   BioMart=.translateChromsForBiomaRt(data$xsome),
                   UCSC=.translateChromsForUCSC(data$xsome),
                   data$xsome)
  ## build GRanges
  ans <- GRanges(seqnames=csomes,
                 ranges=IRanges( ## sign may be reversed
                   start=abs(data$contig_start),
                   end=abs(data$contig_end)),
                 mirna_id = data$mirna_id,
                 strand=data$strand)
  
  ## Filter seqinfo
  .syncSeqlevel(txdb, ans)
}

setGeneric("microRNAs", function(x) standardGeneric("microRNAs"))

## Then set our method
setMethod("microRNAs", "TranscriptDb", function(x){.microRNAs(x)} )



## main function
.tRNAs <- function(txdb){
  require(FDb.UCSC.tRNAs)
  ## get the current package name
  pkgName <- .makePackageName(txdb)
  ## from here we know what the FDB should MUST look like
  fdbName <- sub("TxDb","FDb",pkgName)
  fdbName <- unlist(strsplit(fdbName,"\\."))
  fdbName[5] <- "tRNAs"
  fdbString <- paste(fdbName,collapse=".")
  if(!exists(fdbString)){
    stop("there is no tRNA data available for this organism/source")
  }else{
    ans <- features(eval(parse(text=fdbString)))
  }
  ## Now check active seqs and set the seqlevels
  ans <- .syncSeqlevel(txdb,ans)
  ## now return
  ans
}

setGeneric("tRNAs", function(x) standardGeneric("tRNAs"))

setMethod("tRNAs", "TranscriptDb", function(x){.tRNAs(x)} )

