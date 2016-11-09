### =========================================================================
### The transcripts(), exons(), cds() and promoters() extractors
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_features()
###

### Return a named list of ordinary data frames, 1 per SELECT query.
.extract_features <- function(txdb, proxy_table, mcolumns=character(0),
                              filter=list(), core_columns)
{
    schema_version <- TxDb_schema_version(txdb)
    names(mcolumns) <- TXDB_column2table(mcolumns, from_table=proxy_table,
                                         schema_version=schema_version)
    proxy_column <- orderby <- core_columns[["id"]]

    ## 1st SELECT: extract stuff from the proxy table.
    columns1 <- union(core_columns, mcolumns[names(mcolumns) == proxy_table])
    df1 <- TxDb_SELECT_from_INNER_JOIN(txdb, proxy_table, columns1,
                                       filter=filter, orderby=orderby)

    ## Additional SELECTs: 1 additional SELECT per satellite table with the
    ## exception that the satellite tables that belong to TXDB_SPLICING_BUNDLE
    ## are treated as the virtual single table obtained by LEFT JOIN'ing them
    ## together.
    foreign_columns <- mcolumns[names(mcolumns) != proxy_table]
    bundle_idx <- names(foreign_columns) %in% TXDB_SPLICING_BUNDLE
    names(foreign_columns)[bundle_idx] <- "splicing"
    foreign_columns <- split(foreign_columns, names(foreign_columns))
    satellite_tables <- names(foreign_columns)
    names(satellite_tables) <- satellite_tables
    df_list <- lapply(satellite_tables, function(satellite_table) {
        columns2 <- foreign_columns[[satellite_table]]
        if (length(filter) == 0L) {
            filter2 <- list()
        } else {
            filter2 <- setNames(list(df1[[proxy_column]]), proxy_column)
        }
        if (satellite_table == "splicing") {
            columns2 <- c(proxy_column, columns2)
            if (proxy_table == "transcript") {
                orderby <- c(proxy_column, "exon_rank")
            } else {
                orderby <- c(proxy_column, "_tx_id")
            }
            TxDb_SELECT_from_splicing_bundle(txdb, columns2,
                                             filter=filter2, orderby=orderby)
        } else if (satellite_table == "gene") {
            columns2 <- c(proxy_column, columns2)
            orderby <- c("_tx_id", "gene_id")
            TxDb_SELECT_from_INNER_JOIN(txdb, "gene", columns2,
                                        filter=filter2, orderby=orderby)
        } else {
            stop(satellite_table, ": unsupported satellite table")
        }
    })
    c(setNames(list(df1), proxy_table), df_list)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .extract_features_as_GRanges()
###

.make_DataFrame_from_df_list <- function(df_list)
{
    DF1 <- DataFrame(df_list[[1L]], check.names=FALSE)
    if (length(df_list) == 1L)
        return(DF1)
    proxy_column <- names(DF1)[[1L]]
    ids <- DF1[[proxy_column]]
    DF_list <- lapply(2:length(df_list), function(i) {
        df <- df_list[[i]]
        stopifnot(identical(names(df)[[1L]], proxy_column))
        f <- factor(df[[proxy_column]], levels=ids)
        DataFrame(lapply(df[-1L], function(col) unname(splitAsList(col, f))),
                  check.names=FALSE)
    })
    do.call(DataFrame, c(list(DF1), DF_list, list(check.names=FALSE)))
}

.as_db_columns <- function(columns)
    sub("^(tx_id|exon_id|cds_id)$", "_\\1", columns)

.extract_features_as_GRanges <- function(txdb, proxy_table,
                                         mcolumns=character(0), filter=list())
{
    db_mcolumns <- .as_db_columns(mcolumns)
    names(filter) <- .as_db_columns(names(filter))
    core_columns <- TXDB_table_columns(proxy_table)[TXDB_CORE_TAGS]
    df_list <- .extract_features(txdb, proxy_table, db_mcolumns,
                                 filter, core_columns)
    DF <- .make_DataFrame_from_df_list(df_list)
    ans <- makeGRangesFromDataFrame(DF, seqinfo=get_TxDb_seqinfo0(txdb),
                                    seqnames.field=core_columns[["chrom"]],
                                    start.field=core_columns[["start"]],
                                    end.field=core_columns[["end"]],
                                    strand.field=core_columns[["strand"]])
    mcols(ans) <- setNames(DF[db_mcolumns], mcolumns)
    keep_user_seqlevels_from_TxDb(ans, txdb)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Primary extractors: transcripts(), exons(), cds(), and genes().
###

.dbSchemaHasTxType <- function(txdb){
    txFields <- dbListFields(dbconn(txdb),"transcript")
    "tx_type" %in% txFields
}

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

.extractFromTxDb <- function(txdb, proxy_table,
                             mcolumns=character(0), filter=NULL)
{
    user_mcolumns <- mcolumns
    mcolumns <- translateCols(mcolumns, txdb)
    if (is.null(filter))
        filter <- list()
    names(filter) <- translateCols(names(filter), txdb)
    ans <- .extract_features_as_GRanges(txdb, proxy_table, mcolumns, filter)
    names(mcols(ans)) <- if (is.null(names(user_mcolumns))) user_mcolumns
                         else names(user_mcolumns)
    .assignMetadataList(ans, txdb)
}

setGeneric("transcripts", function(x, ...) standardGeneric("transcripts"))

setMethod("transcripts", "TxDb",
    function(x, columns=c("tx_id", "tx_name"), filter=NULL)
        .extractFromTxDb(x, "transcript", mcolumns=columns, filter=filter)
)

setGeneric("exons", function(x, ...) standardGeneric("exons"))

setMethod("exons", "TxDb",
    function(x, columns="exon_id", filter=NULL)
        .extractFromTxDb(x, "exon", mcolumns=columns, filter=filter)
)

setGeneric("cds", function(x, ...) standardGeneric("cds"))

setMethod("cds", "TxDb",
    function(x, columns="cds_id", filter=NULL)
        .extractFromTxDb(x, "cds", mcolumns=columns, filter=filter)
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
.TxDb.genes <- function(x, columns="gene_id", filter=NULL,
                           single.strand.genes.only=TRUE)
{
    if (!is.character(columns))
        stop("'columns' must be a character vector")
    if (!isTRUEorFALSE(single.strand.genes.only))
        stop("'single.strand.genes.only' must be TRUE or FALSE")
    columns2 <- union(columns, "gene_id")
    tx <- transcripts(x, columns=columns2, filter=filter)

    ## Unroll 'tx' along the 'gene_id' metadata column.
    ## Note that the number of genes per transcript will generally be 1 or 0.
    ## But we also want to handle the situation where it's > 1 which happens
    ## when the same transcript is linked to more than 1 gene (because this
    ## may happen one day and is the reason behind the choice to represent
    ## the 'gene_id' as a CharacterList object instead of a character vector).
    gene_id <- mcols(tx)[ , "gene_id"]
    ngene_per_tx <- elementNROWS(gene_id)
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

    keep_idx <- which(elementNROWS(genes) == 1L)
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

