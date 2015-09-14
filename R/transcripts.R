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
        transcript=c("tx_name", "tx_chrom", "tx_strand", "tx_start", "tx_end",
          "tx_type"),
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
        transcript=c("tx_name", "tx_chrom", "tx_strand", "tx_start", "tx_end",
          "tx_type"),
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


## For converting user arguments FROM the UC style down to what we use
## internally
translateCols <- function(columns, txdb){
    ## preserve any names
    oriColNames <- names(columns)
    ## and save the original column strings
    oriCols <- columns
    
    oriLen <- length(columns) ## not always the same as length(oriCols)
    ## get the available abbreviations as a translation vector (exp)
    names <- .makeColAbbreviations(txdb)
    if(("TXTYPE" %in% columns) && !.dbSchemaHasTxType(txdb)){
        names <- c(names, c(tx_type='TXTYPE'))
    }
    exp <- sub("^_","", names(names))
    names(exp) <- names

    ## Then replace only those IDs that match the UC names
    m <- match(oriCols, names(exp))
    idx = which(!is.na(m))
    columns[idx] <- exp[m[idx]]
    
    if(length(columns) == oriLen && is.null(oriColNames)){
        names(columns) <- oriCols 
    }else if(length(columns) == oriLen && !is.null(oriColNames)){
        names(columns) <- oriColNames         
    }else if(length(columns) != oriLen){
        stop("names were lost in translateCols() helper")
    }    
    columns
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

### In the case of TxDb objects, the distance between 'root_table'
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


.dbSchemaHasTxType <- function(txdb){
    txFields <- dbListFields(dbconn(txdb),"transcript")
    "tx_type" %in% txFields
}

.getSQL_OrderBy <- function(root_table, txdb, what_cols){
    if( ("transcript.tx_type" %in% what_cols) && !.dbSchemaHasTxType(txdb) ){
        ## Only THEN do we adjust the result
        ## Then drop that string
        what_cols <- what_cols[!(what_cols %in% "transcript.tx_type")]
        what <- paste(what_cols, collapse=", ")
    }else{
        what <- paste(what_cols, collapse=", ")
    }
    what
}

.getSQL_What <- function(root_table, txdb, what_cols){
    if( ("transcript.tx_type" %in% what_cols) && !.dbSchemaHasTxType(txdb) ){
        ## Only THEN do we adjust the result
        what <- what_cols    
        old <- "transcript.tx_type"  ###paste0(root_table,'.tx_type')
        new <- "NULL AS tx_type"
        ## then replace old with new
        what[what %in% old] <- new
        what <- paste(what, collapse=", ")
    }else{
        what <- paste(what_cols, collapse=", ")
    }
    what
}

.extractData <- function(root_table, txdb, what_cols, child_tables, vals,
                         orderby_cols=NULL)
{
    ## SQL_what <- paste(.asQualifiedColnames(root_table, what_cols),
    ##                   collapse=", ")
    what <- .asQualifiedColnames(root_table, what_cols)
    SQL_what <- .getSQL_What(root_table, txdb, what)
    SQL_from <- .joinRootToChildTables(root_table, child_tables)
    SQL_where <- .sqlWhereIn(vals)
    if (length(orderby_cols) == 0L) {
        SQL_orderby <- ""
    } else {
        orderby <- .asQualifiedColnames(root_table, orderby_cols)
        SQL_orderby <- paste("ORDER BY", paste(orderby, collapse=", "))
    }
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
    ## Ordering by '_tx_id' guarantees that the objects returned by
    ## transcripts() and exonsBy() are parallel.
    #orderby_cols <- CORECOLS[c("chrom", "strand", "start", "end")]
    orderby_cols <- CORECOLS["id"]
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
    what <- c(primary_key,
                   .asQualifiedColnames(root_table, child_columns))
    
##    SQL_what <- paste(what, collapse=", ")
    SQL_what <- .getSQL_What(root_table, txdb, what)
    SQL_from <- .joinPrimaryKeyToChildTable(root_table, child_table)
    vals <- list(ids)
    names(vals) <- primary_key
    SQL_where <- .sqlWhereIn(vals)
    ## SQL_orderby <- SQL_what ## similar fix here to the SQL_what
    SQL_orderby <- .getSQL_OrderBy(root_table, txdb, what)
    SQL <- paste("SELECT DISTINCT", SQL_what, "FROM", SQL_from,
                 SQL_where, "ORDER BY", SQL_orderby)
    ans <- dbEasyQuery(dbconn(txdb), SQL)
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

## me a named list from the metadata data.frame
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

## helper to translate back to what is expected from seqinfo()
.translateToSeqInfo <- function(txdb, x){
    tr <- load_chrominfo(txdb, set.col.class=TRUE)$chrom[txdb$new2old]
    names(tr) <- txdb$.chrom    
    idx <- match(x, tr)
    names(tr)[idx]
}

.extractFeatureRowsAsGRanges <- function(root_table, txdb, vals, columns)
{
    ## 1st translate columns from UC format to LC format
    columns <- translateCols(columns, txdb)
    ## Then proceed with checking
    CORECOLS <- .DBDESC[[root_table]]$CORECOLS
    assigned_columns <- .assignColToClosestTable(root_table, columns) 
    root_columns <- assigned_columns[[root_table]]
    ## Extract the data from the db.
    root_data <- .extractRootData(root_table, txdb, vals, root_columns)
    ## seqnames may be out of sync with expected results.  Massage back.
    root_data[[CORECOLS["chrom"]]] <- .translateToSeqInfo(txdb, 
                                          root_data[[CORECOLS["chrom"]]])
    child_data <- .extractChildData(root_table, txdb,
                          root_data[[CORECOLS["id"]]], assigned_columns)
    ## Construct the GRanges object and return it.
    ans_seqinfo <- seqinfo(txdb)
    ans_seqnames <- factor(root_data[[CORECOLS["chrom"]]],
                           levels=seqlevels(ans_seqinfo))
    ans_ranges <- IRanges(start=root_data[[CORECOLS["start"]]],
                          end=root_data[[CORECOLS["end"]]])
    ans_strand <- strand(root_data[[CORECOLS["strand"]]])

    activeNames <- names(.isActiveSeq(txdb))[.isActiveSeq(txdb)]
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



## this helper is just to get the .isActiveSeq vector, but to have it
## named based on the original database seqnames...
## This is important for places where .isActiveSeq() needs to be used
## as part of a database query instead of as part of a external
## representation.
.baseNamedActiveSeqs <- function(txdb){
    trueNames <- load_chrominfo(txdb, set.col.class=TRUE)$chrom
    actSqs <- .isActiveSeq(txdb)
    names(actSqs) <- trueNames[txdb$new2old] ## limit result to these.
    actSqs
}

## This is used to create a list from the .isActiveSeq slot
## !!!
## TODO: I think that this helper is screwing up the vals values in the methods
.makeActiveSeqsList <- function(type, txdb){
    actSqs <- .baseNamedActiveSeqs(txdb)
    keepSeqs <- names(actSqs)[actSqs]
    res <- list(keepSeqs)
    names(res) <- type
    res
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Primary extractors: transcripts(), exons(), cds(), and genes().
###
### TODO: Rename the 'vals' arg -> 'filter' so it's consistent with the
### 'filters' arg of makeTxDbFromBiomart() which we should also rename
### 'filter'.
### Other proposal: rename the 'columns' arg -> 'colnames'
###

setGeneric("transcripts", function(x, ...) standardGeneric("transcripts"))

## TODOS: change defaults (WILL break many examples!)
## TODO: change defaults from c("tx_id", "tx_name") to: c("TXID", "TXNAME") 
setMethod("transcripts", "TxDb",
    function(x, vals=NULL, columns=c("tx_id", "tx_name")){
        vals = c(vals, .makeActiveSeqsList("tx_chrom", x))
        .extractFeatureRowsAsGRanges("transcript", x, vals, columns)
      }
)

setGeneric("exons", function(x, ...) standardGeneric("exons"))

## TODO: change defaults from c("exon_id") to: c("EXONID") 
setMethod("exons", "TxDb",
    function(x, vals=NULL, columns="exon_id"){
        vals = c(vals, .makeActiveSeqsList("exon_chrom", x))
        .extractFeatureRowsAsGRanges("exon", x, vals, columns)
        }
)

setGeneric("cds", function(x, ...) standardGeneric("cds"))

## TODO: change defaults from c("cds_id") to: c("CDSID") 
setMethod("cds", "TxDb",
    function(x, vals=NULL, columns="cds_id"){
        vals = c(vals, .makeActiveSeqsList("cds_chrom", x))
        .extractFeatureRowsAsGRanges("cds", x, vals, columns)
        }
)

setGeneric("genes", function(x, ...) standardGeneric("genes"))

.relist_col <- function(x, skeleton)
{
   if (is.list(x) || (is(x, "List") && !is(x, "Ranges")))
       return(IRanges:::regroupBySupergroup(x, skeleton))
   relist(x, skeleton)
}

### Used in SplicingGraphs! Move this to IRanges (and check similarity with
### "collapse" method for DataFrame objects).
.collapse_df <- function(df, skeleton)
{
    ## FIXME: endoapply() on a DataFrame object is broken when applying
    ## a function 'FUN' that modifies the nb of rows. Furthermore, the
    ## returned object passes validation despite being broken! Fix it
    ## in IRanges.
    ans <- endoapply(df, function(x) unique(.relist_col(x, skeleton)))
    ## Fix the broken DataFrame returned by endoapply().
    ans@nrows <- length(skeleton)
    ans@rownames <- NULL
    ans
}

### If 'single.strand.genes.only' is TRUE (the default), then genes that
### have exons located on both strands of the same chromosome, or on 2
### different chromosomes are dropped. In that case, the genes are returned
### in a GRanges object. Otherwise, they're returned in a GRangesList object
### with the metadata columns requested thru 'columns' set at the top level.
.TxDb.genes <- function(x, vals=NULL, columns="gene_id",
                        single.strand.genes.only=TRUE)
{
    if (!is.character(columns))
        stop("'columns' must be a character vector")
    if (!isTRUEorFALSE(single.strand.genes.only))
        stop("'single.strand.genes.only' must be TRUE or FALSE")
    columns2 <- union(columns, "gene_id")
    tx <- transcripts(x, vals=vals, columns=columns2)

    ## Unroll 'tx' along the 'gene_id' metadata column.
    ## Note that the number of genes per transcript will generally be 1 or 0.
    ## But we also want to handle the situation where it's > 1 which happens
    ## when the same transcript is linked to more than 1 gene (because this
    ## may happen one day and is the reason behind the choice to represent
    ## the 'gene_id' as a CharacterList object instead of a character vector).
    gene_id <- mcols(tx)[ , "gene_id"]
    ngene_per_tx <- elementLengths(gene_id)
    tx <- tx[rep.int(seq_along(ngene_per_tx), ngene_per_tx)]
    mcols(tx)$gene_id <- unlist(gene_id, use.names=FALSE)

    ## Split 'tx' by gene.
    tx_by_gene <- split(tx, mcols(tx)$gene_id)

    ## Turn inner mcols into outter mcols by relist()'ing them.
    inner_mcols <- mcols(tx_by_gene@unlistData)[columns]
    mcols(tx_by_gene@unlistData) <- NULL
    outter_mcols <- .collapse_df(inner_mcols, tx_by_gene)
    gene_id <- outter_mcols$gene_id
    if (!is.null(gene_id)) {
        stopifnot(identical(names(tx_by_gene), as.character(gene_id)))
        outter_mcols$gene_id <- names(tx_by_gene)
    }
    mcols(tx_by_gene) <- outter_mcols

    ## Compute the gene ranges.
    genes <- range(tx_by_gene)

    if (!single.strand.genes.only)
        return(genes)

    keep_idx <- which(elementLengths(genes) == 1L)
    genes <- genes[keep_idx]
    ans <- unlist(genes, use.names=FALSE)
    mcols(ans) <- mcols(genes)
    names(ans) <- names(genes)
    ans
}

setMethod("genes", "TxDb", .TxDb.genes)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "promoters" method
###
### generic is in IRanges
###

setMethod("promoters", "TxDb",
    function(x, upstream=2000, downstream=200, ...)
    {
        gr <- transcripts(x, ...)
        promoters(gr, upstream=upstream, downstream=downstream)
    }
) 

