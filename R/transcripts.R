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
                                         mcolumns=character(0), filter=list(),
                                         use.names=FALSE)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    db_mcolumns <- db_mcolumns0 <- .as_db_columns(mcolumns)
    proxy_columns <- TXDB_table_columns(proxy_table)
    if (use.names) {
        if ("name" %in% names(proxy_columns)) {
            proxy_name_column <- proxy_columns[["name"]]
            if (!(proxy_name_column %in% db_mcolumns0))
                db_mcolumns <- c(db_mcolumns0, proxy_name_column)
        } else {
            warning(wmsg("no column in '", proxy_table, "' table ",
                         "to retrieve the feature names from"))
            use.names <- FALSE
        }
    }
    names(filter) <- .as_db_columns(names(filter))
    core_columns <- proxy_columns[TXDB_CORE_COLTAGS]
    df_list <- .extract_features(txdb, proxy_table, db_mcolumns,
                                 filter, core_columns)
    DF <- .make_DataFrame_from_df_list(df_list)
    ans <- makeGRangesFromDataFrame(DF, seqinfo=get_TxDb_seqinfo0(txdb),
                                    seqnames.field=core_columns[["chrom"]],
                                    start.field=core_columns[["start"]],
                                    end.field=core_columns[["end"]],
                                    strand.field=core_columns[["strand"]])
    if (use.names) {
        ans_names <- DF[ , proxy_name_column]
        ans_names[is.na(ans_names)] <- ""  # replace NAs with empty strings
        names(ans) <- ans_names
    }
    mcols(ans) <- setNames(DF[db_mcolumns0], mcolumns)
    set_user_seqlevels_and_genome(ans, txdb)
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
                             mcolumns=character(0), filter=NULL,
                             use.names=FALSE)
{
    user_mcolumns <- mcolumns
    mcolumns <- translateCols(mcolumns, txdb)
    if (is.null(filter))
        filter <- list()
    names(filter) <- translateCols(names(filter), txdb)
    ans <- .extract_features_as_GRanges(txdb, proxy_table, mcolumns, filter,
                                        use.names)
    names(mcols(ans)) <- if (is.null(names(user_mcolumns))) user_mcolumns
                         else names(user_mcolumns)
    .assignMetadataList(ans, txdb)
}

#' Extract genomic features from a TxDb-like object
#'
#' Generic functions to extract genomic features from a TxDb-like object.  This
#' page documents the methods for \link{TxDb} objects only.
#'
#' These are the main functions for extracting transcript information from a
#' \link{TxDb}-like object. These methods can restrict the output based on
#' categorical information. To restrict the output based on interval
#' information, use the \code{\link{transcriptsByOverlaps}},
#' \code{\link{exonsByOverlaps}}, and \code{\link{cdsByOverlaps}} functions.
#'
#' The \code{promoters} function computes user-defined promoter regions for the
#' transcripts in a \link{TxDb}-like object. The return object is a
#' \code{GRanges} of promoter regions around the transcription start site the
#' span of which is defined by \code{upstream} and \code{downstream}.  For
#' additional details on how the promoter range is computed and the handling of
#' \code{+} and \code{-} strands see \code{?`promoters,GRanges-method`}.
#'
#' @aliases transcripts transcripts,TxDb-method exons exons,TxDb-method cds
#' cds,TxDb-method genes genes,TxDb-method promoters promoters,TxDb-method
#' @param x A \link{TxDb} object.
#' @param ...  For the \code{transcripts}, \code{exons}, \code{cds}, and
#' \code{genes} generic functions: arguments to be passed to methods.
#'
#' For the \code{promoters} method for \link{TxDb} objects: arguments to be
#' passed to the internal call to \code{transcripts}.
#' @param columns Columns to include in the output.  Must be \code{NULL} or a
#' character vector as given by the \code{columns} method. With the following
#' restrictions: \itemize{ \item \code{"TXCHROM"} and \code{"TXSTRAND"} are not
#' allowed for \code{transcripts}.  \item \code{"EXONCHROM"} and
#' \code{"EXONSTRAND"} are not allowed for \code{exons}.  \item
#' \code{"CDSCHROM"} and \code{"CDSSTRAND"} are not allowed for \code{cds}.  }
#' If the vector is named, those names are used for the corresponding column in
#' the element metadata of the returned object.
#' @param filter Either \code{NULL} or a named list of vectors to be used to
#' restrict the output. Valid names for this list are: \code{"gene_id"},
#' \code{"tx_id"}, \code{"tx_name"}, \code{"tx_chrom"}, \code{"tx_strand"},
#' \code{"exon_id"}, \code{"exon_name"}, \code{"exon_chrom"},
#' \code{"exon_strand"}, \code{"cds_id"}, \code{"cds_name"},
#' \code{"cds_chrom"}, \code{"cds_strand"} and \code{"exon_rank"}.
#' @param use.names \code{TRUE} or \code{FALSE}. If \code{TRUE}, the feature
#' names are set as the names of the returned object, with NAs being replaced
#' with empty strings.
#' @param single.strand.genes.only \code{TRUE} or \code{FALSE}.  If \code{TRUE}
#' (the default), then genes are returned in a \link[GenomicRanges]{GRanges}
#' object and those genes that cannot be represented by a single genomic range
#' (because they have exons located on both strands of the same reference
#' sequence or on more than one reference sequence) are dropped with a message.
#'
#' If \code{FALSE}, then all the genes are returned in a
#' \link[GenomicRanges]{GRangesList} object with the columns specified thru the
#' \code{columns} argument set as \emph{top level} metadata columns. (Please
#' keep in mind that the \emph{top level} metadata columns of a
#' \link[GenomicRanges]{GRangesList} object are not displayed by the
#' \code{show()} method.)
#' @param upstream For \code{promoters} : An \code{integer(1)} value indicating
#' the number of bases upstream from the transcription start site. For
#' additional details see \code{?`promoters,GRanges-method`}.
#' @param downstream For \code{promoters} : An \code{integer(1)} value
#' indicating the number of bases downstream from the transcription start site.
#' For additional details see \code{?`promoters,GRanges-method`}.
#' @return A \link[GenomicRanges]{GRanges} object. The only exception being
#' when \code{genes} is used with \code{single.strand.genes.only=FALSE}, in
#' which case a \link[GenomicRanges]{GRangesList} object is returned.
#' @author M. Carlson, P. Aboyoun and H. Pagès
#' @seealso \itemize{ \item \code{\link{transcriptsBy}} and
#' \code{\link{transcriptsByOverlaps}} for more ways to extract genomic
#' features from a \link{TxDb}-like object.
#'
#' \item \code{\link{transcriptLengths}} for extracting the transcript lengths
#' (and other metrics) from a \link{TxDb} object.
#'
#' \item \code{\link{exonicParts}} and \code{\link{intronicParts}} for
#' extracting non-overlapping exonic or intronic parts from a TxDb-like object.
#'
#' \item \code{\link{extendExonsIntoIntrons}} for extending exons into their
#' adjacent introns.
#'
#' \item \code{\link{extractTranscriptSeqs}} for extracting transcript (or CDS)
#' sequences from reference sequences.
#'
#' \item \code{\link{coverageByTranscript}} for computing coverage by
#' transcript (or CDS) of a set of ranges.
#'
#' \item \link[GenomicFeatures]{select-methods} for how to use the simple
#' "select" interface to extract information from a \link{TxDb} object.
#'
#' \item \code{\link{microRNAs}} and \code{\link{tRNAs}} for extracting
#' microRNA or tRNA genomic ranges from a \link{TxDb} object.
#'
#' \item \code{\link{id2name}} for mapping \link{TxDb} internal ids to external
#' names for a given feature type.
#'
#' \item The \link{TxDb} class.  }
#' @keywords methods
#' @examples
#'
#' txdb_file <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                          package="GenomicFeatures")
#' txdb <- loadDb(txdb_file)
#'
#' ## ---------------------------------------------------------------------
#' ## transcripts()
#' ## ---------------------------------------------------------------------
#'
#' tx1 <- transcripts(txdb)
#' tx1
#'
#' transcripts(txdb, use.names=TRUE)
#' transcripts(txdb, columns=NULL, use.names=TRUE)
#'
#' filter <- list(tx_chrom = c("chr3", "chr5"), tx_strand = "+")
#' tx2 <- transcripts(txdb, filter=filter)
#' tx2
#'
#' ## Sanity checks:
#' stopifnot(
#'   identical(mcols(tx1)$tx_id, seq_along(tx1)),
#'   identical(tx2, tx1[seqnames(tx1) == "chr3" & strand(tx1) == "+"])
#' )
#'
#' ## ---------------------------------------------------------------------
#' ## exons()
#' ## ---------------------------------------------------------------------
#'
#' exons(txdb, columns=c("EXONID", "TXNAME"),
#'             filter=list(exon_id=1))
#' exons(txdb, columns=c("EXONID", "TXNAME"),
#'             filter=list(tx_name="uc009vip.1"))
#'
#' ## ---------------------------------------------------------------------
#' ## genes()
#' ## ---------------------------------------------------------------------
#'
#' genes(txdb)  # a GRanges object
#' cols <- c("tx_id", "tx_chrom", "tx_strand",
#'           "exon_id", "exon_chrom", "exon_strand")
#' ## By default, genes are returned in a GRanges object and those that
#' ## cannot be represented by a single genomic range (because they have
#' ## exons located on both strands of the same reference sequence or on
#' ## more than one reference sequence) are dropped with a message:
#' single_strand_genes <- genes(txdb, columns=cols)
#'
#' ## Because we've returned single strand genes only, the "tx_chrom"
#' ## and "exon_chrom" metadata columns are guaranteed to match
#' ## 'seqnames(single_strand_genes)':
#' stopifnot(identical(as.character(seqnames(single_strand_genes)),
#'                     as.character(mcols(single_strand_genes)$tx_chrom)))
#' stopifnot(identical(as.character(seqnames(single_strand_genes)),
#'                     as.character(mcols(single_strand_genes)$exon_chrom)))
#'
#' ## and also the "tx_strand" and "exon_strand" metadata columns are
#' ## guaranteed to match 'strand(single_strand_genes)':
#' stopifnot(identical(as.character(strand(single_strand_genes)),
#'                     as.character(mcols(single_strand_genes)$tx_strand)))
#' stopifnot(identical(as.character(strand(single_strand_genes)),
#'                     as.character(mcols(single_strand_genes)$exon_strand)))
#'
#' all_genes <- genes(txdb, columns=cols, single.strand.genes.only=FALSE)
#' all_genes  # a GRangesList object
#' multiple_strand_genes <- all_genes[elementNROWS(all_genes) >= 2]
#' multiple_strand_genes
#' mcols(multiple_strand_genes)
#'
#' ## ---------------------------------------------------------------------
#' ## promoters()
#' ## ---------------------------------------------------------------------
#'
#' ## This:
#' promoters(txdb, upstream=100, downstream=50)
#' ## is equivalent to:
#' promoters(transcripts(txdb, use.names=TRUE), upstream=100, downstream=50)
#'
#' ## Extra arguments are passed to transcripts(). So this:
#' columns <- c("tx_name", "gene_id")
#' promoters(txdb, upstream=100, downstream=50, columns=columns)
#' ## is equivalent to:
#' promoters(transcripts(txdb, columns=columns, use.names=TRUE),
#'           upstream=100, downstream=50)
#'
#' @export transcripts
setGeneric("transcripts", function(x, ...) standardGeneric("transcripts"))

setMethod("transcripts", "TxDb",
    function(x, columns=c("tx_id", "tx_name"), filter=NULL, use.names=FALSE)
        .extractFromTxDb(x, "transcript", mcolumns=columns, filter=filter,
                                          use.names=use.names)
)

setGeneric("exons", function(x, ...) standardGeneric("exons"))

setMethod("exons", "TxDb",
    function(x, columns="exon_id", filter=NULL, use.names=FALSE)
        .extractFromTxDb(x, "exon", mcolumns=columns, filter=filter,
                                    use.names=use.names)
)

setGeneric("cds", function(x, ...) standardGeneric("cds"))

setMethod("cds", "TxDb",
    function(x, columns="cds_id", filter=NULL, use.names=FALSE)
        .extractFromTxDb(x, "cds", mcolumns=columns, filter=filter,
                                   use.names=use.names)
)

setGeneric("genes", function(x, ...) standardGeneric("genes"))

.relist_col <- function(x, skeleton)
{
   if (is.list(x) || (is(x, "List") && !is(x, "IntegerRanges")))
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

### If 'single.strand.genes.only' is TRUE (the default), then genes are
### returned in a GRanges object and those genes that cannot be represented
### by a single genomic range (because they have exons located on both strands
### of the same reference sequence or on more than one reference sequence)
### are dropped with a message. Otherwise, they're returned in a GRangesList
### object with the metadata columns requested thru 'columns' set at the top
### level.
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

    is_single_range_gene <- elementNROWS(genes) == 1L
    nb_multi_range_genes <- sum(!is_single_range_gene)
    if (nb_multi_range_genes != 0L) {
        if (nb_multi_range_genes == 1L) {
            what <- "gene was dropped because it has"
        } else {
            what <- "genes were dropped because they have"
        }
        message("  ", wmsg(nb_multi_range_genes, " ", what, " exons located ",
                "on both strands of the same reference sequence or on more ",
                "than one reference sequence, so cannot be represented by a ",
                "single genomic range."), "\n  ",
                wmsg("Use 'single.strand.genes.only=FALSE' to get all the ",
                "genes in a GRangesList object, or use suppressMessages() ",
                "to suppress this message."))
    }
    keep_idx <- which(is_single_range_gene)
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
    function(x, upstream=2000, downstream=200, use.names=TRUE, ...)
    {
        gr <- transcripts(x, ..., use.names=use.names)
        promoters(gr, upstream=upstream, downstream=downstream)
    }
)

