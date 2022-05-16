### =========================================================================
### The "features grouped by" extractors
### -------------------------------------------------------------------------


.set_group_names <- function(ans, use.names, txdb, by)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (!use.names)
        return(ans)
    if (by == "gene")
        stop("'use.names=TRUE' cannot be used when grouping by gene")
    names(ans) <- id2name(txdb, feature.type=by)[names(ans)]
    if (any(is.na(names(ans))) || any(duplicated(names(ans))))
        warning("some group names are NAs or duplicated")
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### transcriptsBy(), exonsBy() and cdsBy().
###

#' Extract and group genomic features of a given type from a TxDb-like object
#'
#' Generic functions to extract genomic features of a given type grouped based
#' on another type of genomic feature.  This page documents the methods for
#' \link{TxDb} objects only.
#'
#' These functions return a \link[GenomicRanges]{GRangesList} object where the
#' ranges within each of the elements are ordered according to the following
#' rule:
#'
#' When using \code{exonsBy} or \code{cdsBy} with \code{by = "tx"}, the
#' returned exons or CDS are ordered by ascending rank for each transcript,
#' that is, by their position in the transcript.  In all other cases, the
#' ranges will be ordered by chromosome, strand, start, and end values.
#'
#' @aliases transcriptsBy transcriptsBy,TxDb-method exonsBy exonsBy,TxDb-method
#' cdsBy cdsBy,TxDb-method intronsByTranscript intronsByTranscript,TxDb-method
#' fiveUTRsByTranscript fiveUTRsByTranscript,TxDb-method threeUTRsByTranscript
#' threeUTRsByTranscript,TxDb-method
#' @param x A \link{TxDb} object.
#' @param ... Arguments to be passed to or from methods.
#' @param by One of \code{"gene"}, \code{"exon"}, \code{"cds"} or \code{"tx"}.
#' Determines the grouping.
#' @param use.names Controls how to set the names of the returned
#' \link[GenomicRanges]{GRangesList} object.  These functions return all the
#' features of a given type (e.g.  all the exons) grouped by another feature
#' type (e.g. grouped by transcript) in a \link[GenomicRanges]{GRangesList}
#' object.  By default (i.e. if \code{use.names} is \code{FALSE}), the names of
#' this \link[GenomicRanges]{GRangesList} object (aka the group names) are the
#' internal ids of the features used for grouping (aka the grouping features),
#' which are guaranteed to be unique.  If \code{use.names} is \code{TRUE}, then
#' the names of the grouping features are used instead of their internal ids.
#' For example, when grouping by transcript (\code{by="tx"}), the default group
#' names are the transcript internal ids (\code{"tx_id"}). But, if
#' \code{use.names=TRUE}, the group names are the transcript names
#' (\code{"tx_name"}).  Note that, unlike the feature ids, the feature names
#' are not guaranteed to be unique or even defined (they could be all
#' \code{NA}s). A warning is issued when this happens.  See
#' \code{?\link{id2name}} for more information about feature internal ids and
#' feature external names and how to map the formers to the latters.
#'
#' Finally, \code{use.names=TRUE} cannot be used when grouping by gene
#' \code{by="gene"}. This is because, unlike for the other features, the gene
#' ids are external ids (e.g. Entrez Gene or Ensembl ids) so the db doesn't
#' have a \code{"gene_name"} column for storing alternate gene names.
#' @return A \link[GenomicRanges]{GRangesList} object.
#' @author M. Carlson, P. Aboyoun and H. Pagès
#' @seealso \itemize{ \item \code{\link{transcripts}} and
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
#' sequences from chromosome sequences.
#'
#' \item \code{\link{coverageByTranscript}} for computing coverage by
#' transcript (or CDS) of a set of ranges.
#'
#' \item \link[GenomicFeatures]{select-methods} for how to use the simple
#' "select" interface to extract information from a \link{TxDb} object.
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
#' ## Get the transcripts grouped by gene:
#' transcriptsBy(txdb, "gene")
#'
#' ## Get the exons grouped by gene:
#' exonsBy(txdb, "gene")
#'
#' ## Get the CDS grouped by transcript:
#' cds_by_tx0 <- cdsBy(txdb, "tx")
#' ## With more informative group names:
#' cds_by_tx1 <- cdsBy(txdb, "tx", use.names=TRUE)
#' ## Note that 'cds_by_tx1' can also be obtained with:
#' names(cds_by_tx0) <- id2name(txdb, feature.type="tx")[names(cds_by_tx0)]
#' stopifnot(identical(cds_by_tx0, cds_by_tx1))
#'
#' ## Get the introns grouped by transcript:
#' intronsByTranscript(txdb)
#'
#' ## Get the 5' UTRs grouped by transcript:
#' fiveUTRsByTranscript(txdb)
#' fiveUTRsByTranscript(txdb, use.names=TRUE)  # more informative group names
#'
NULL

.join_genes_and_Rdf <- function(genes, Rdf, using)
{
    gene_id <- genes[ , "gene_id"]
    gene_id <- factor(gene_id, levels=unique(gene_id))
    join_idx <- match(genes[ , using], Rdf[ , using])
    Rdf <- S4Vectors:::extract_data_frame_rows(Rdf, join_idx)
    cbind(gene_id=gene_id, Rdf)
}

### 'columns' must be a named vector of db columns where the names are user
### columns.
.add_prefix_to_user_columns <- function(columns, proxy_table)
{
    if (proxy_table == "transcript") {
        prefix <- "tx"
    } else {
        prefix <- proxy_table
    }
    user_columns <- names(columns)
    has_noname <- user_columns %in% c(NA_character_, "")
    prefix_idx <- which(!has_noname)
    user_columns[prefix_idx] <- paste(prefix, user_columns[prefix_idx], sep="_")
    user_columns[has_noname] <- columns[has_noname]
    names(columns) <- user_columns
    columns
}

### 'columns' must be a named character vector with the names being used to
### rename the columns.
.split_df_into_GRL <- function(txdb, df, columns, by, use.names)
{
    rownames(df) <- NULL
    f <- df[[1L]]
    gr <- makeGRangesFromDataFrame(df[-1L], seqinfo=get_TxDb_seqinfo0(txdb))
    gr_mcols <- setNames(df[columns], names(columns))
    mcols(gr) <- DataFrame(gr_mcols, check.names=FALSE)
    grl <- split(gr, f)
    grl <- set_user_seqlevels_and_genome(grl, txdb)
    grl <- .set_group_names(grl, use.names, txdb, by)
    .assignMetadataList(grl, txdb)
}

.extract_features_by_gene <- function(txdb, proxy_table, columns)
{
    join_column <- columns[["id"]]
    orderby <- c("gene_id", join_column)

    ## We use 2 SELECT queries and join the results ourselves. This join
    ## is actually faster than the SQL JOIN.

    ## 1st SELECT query: get the 2-column 'genes' data frame.
    if (proxy_table == "transcript") {
        table <- "gene"
    } else {
        table <- "splicing"
    }
    genes <- TxDb_SELECT_from_INNER_JOIN(txdb, table, orderby,
                                         orderby=orderby)
    if (proxy_table == "cds") {
        ## Proxy column is from the "splicing" table, not from the "cds"
        ## table (which was not even involved in the JOIN in the first
        ## place), so it can contain NAs. Remove these rows.
        keep_me <- !is.na(genes[[2L]])
        genes <- S4Vectors:::extract_data_frame_rows(genes, keep_me)
    }

    ## 2nd SELECT query.
    Rdf <- TxDb_SELECT_from_INNER_JOIN(txdb, proxy_table, columns)

    ## Join the results.
    .join_genes_and_Rdf(genes, Rdf, join_column)
}

.extract_features_by <- function(txdb, proxy_table, by, use.names=FALSE)
{
    tags <- c(TXDB_CORE_COLTAGS, "name")
    columns <- TXDB_table_columns(proxy_table)[tags]
    if (by == "gene") {
        df <- .extract_features_by_gene(txdb, proxy_table, columns)
        ## Inner metadata columns of the returned GRangesList object.
        mcolumns <- columns[c("id", "name")]
    } else {
        if (proxy_table == "transcript") {
            ## Extract transcripts grouped by exon or CDS.
            by_column <- TXDB_table_columns(by)[["id"]]
            orderby <- c(by_column, columns[["id"]], "exon_rank")
        } else {
            ## Extract exons or CDS grouped by transcript.
            by_column <- TXDB_table_columns("transcript")[["id"]]
            orderby <- c(by_column, "exon_rank")
        }
        columns <- c(by_column, "exon_rank", columns)
        df <- TxDb_SELECT_from_splicing_bundle(txdb, columns,
                                               orderby=orderby,
                                               cds_join_type="INNER")
        ## Inner metadata columns of the returned GRangesList object.
        mcolumns <- c(columns[c("id", "name")], "exon_rank")
    }
    mcolumns <- .add_prefix_to_user_columns(mcolumns, proxy_table)
    .split_df_into_GRL(txdb, df, mcolumns, by, use.names)
}

setGeneric("transcriptsBy", signature="x",
    function(x, by=c("gene", "exon", "cds"), ...)
        standardGeneric("transcriptsBy")
)

setMethod("transcriptsBy", "TxDb",
    function(x, by=c("gene", "exon", "cds"), use.names=FALSE)
    {
        by <- match.arg(by)
        .extract_features_by(x, "transcript", by, use.names=use.names)
    }
)

setGeneric("exonsBy", signature="x",
    function(x, by=c("tx", "gene"), ...) standardGeneric("exonsBy")
)

setMethod("exonsBy", "TxDb",
    function(x, by=c("tx", "gene"), use.names=FALSE)
    {
        by <- match.arg(by)
        .extract_features_by(x, "exon", by, use.names=use.names)
    }
)

setGeneric("cdsBy", signature="x",
    function(x, by=c("tx", "gene"), ...) standardGeneric("cdsBy")
)

setMethod("cdsBy", "TxDb",
    function(x, by=c("tx", "gene"), use.names=FALSE)
    {
        by <- match.arg(by)
        .extract_features_by(x, "cds", by, use.names=use.names)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### intronsByTranscript()
###

setGeneric("intronsByTranscript",
    function(x, ...) standardGeneric("intronsByTranscript")
)

setMethod("intronsByTranscript", "TxDb",
    function(x, use.names=FALSE)
    {
        tx <- transcripts(x)
        exn <- exonsBy(x)
        tx <- tx[match(names(exn), mcols(tx)[,"tx_id"])]
        ans <- psetdiff(tx, exn)
        .set_group_names(ans, use.names, x, "tx")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fiveUTRsByTranscript() and threeUTRsByTranscript()
###

.getSplicingsForTranscriptsWithCDSs <- function(txdb)
{
    ans <- load_splicings(txdb)
    ids <- unique(ans$tx_id[!is.na(ans$cds_id)])
    S4Vectors:::extract_data_frame_rows(ans, ans$tx_id %in% ids)
}

### 'tx_id': character or integer vector with runs of identical elements (one
### run per transcript, and, within each run, one element per exon).
### 'exons_with_cds': integer vector containing the indices of the elements
### in 'tx_id' corresponding to exons with a CDS. The indices must be sorted
### in ascending order.
### Returns for each transcript the indices of the first exon with a CDS, plus
### the indices of the preceding exons in the transcript.
### Note that, for each transcript, exons that have a 5' UTR are all the exons
### before the first exon with a CDS, including the first exon with a CDS (even
### though this one might actually have a 0-width 5' UTR but we will take care
### of this later).
.exons_with_5utr <- function(tx_id, exons_with_cds)
{
    tx_id_rle <- Rle(tx_id)
    nexon <- runLength(tx_id_rle)
    ntx <- length(nexon)
    pseudo_tx_id <- rep.int(seq_len(ntx), nexon)
    first_exon_with_cds <- integer(ntx)
    exons_with_cds <- rev(exons_with_cds)
    first_exon_with_cds[pseudo_tx_id[exons_with_cds]] <- exons_with_cds
    offset <- cumsum(c(0L, nexon[-length(nexon)]))
    lengths <- first_exon_with_cds - offset
    sequence(lengths, from=offset+1L)
}

### 'tx_id', 'exons_with_cds': same as for .exons_with_5utr().
### Returns for each transcript the indices of the last exon with a CDS, plus
### the indices of the following exons in the transcript.
### Note that, for each transcript, exons that have a 3' UTR are all the exons
### after the last exon with a CDS, including the last exon with a CDS (even
### though this one might actually have a 0-width 3' UTR but we will take care
### of this later).
.exons_with_3utr <- function(tx_id, exons_with_cds)
{
    tx_id_rle <- Rle(tx_id)
    nexon <- runLength(tx_id_rle)
    ntx <- length(nexon)
    pseudo_tx_id <- rep.int(seq_len(ntx), nexon)
    last_exon_with_cds <- integer(ntx)
    last_exon_with_cds[pseudo_tx_id[exons_with_cds]] <- exons_with_cds
    lengths <- cumsum(nexon) - last_exon_with_cds + 1L
    sequence(lengths, from=last_exon_with_cds)
}

.makeUTRsByTranscript <- function(txdb, splicings, utr_start, utr_end)
{
    seqinfo0 <- get_TxDb_seqinfo0(txdb)
    cols <- paste0("exon_", c("id", "name", "rank"))
    gr <- GRanges(seqnames = factor(
                     splicings$exon_chrom,
                     levels =seqlevels(seqinfo0)),
                  ranges = IRanges(start = utr_start, end = utr_end),
                  strand = strand(splicings$exon_strand),
                  splicings[cols],
                  seqinfo = seqinfo0)
    idx <- width(gr) != 0L  # drop 0-width UTRs
    grl <- split(gr[idx], splicings$tx_id[idx])
    set_user_seqlevels_and_genome(grl, txdb)
}

.make5UTRsByTranscript <- function(txdb, splicings, use.names=FALSE)
{
    exons_with_cds <- which(!is.na(splicings$cds_id))
    idx <- .exons_with_5utr(splicings$tx_id, exons_with_cds)
    splicings <- S4Vectors:::extract_data_frame_rows(splicings, idx)

    ## Compute the UTR starts/ends.
    utr_start <- splicings$exon_start
    utr_end <- splicings$exon_end
    idx1 <- !is.na(splicings$cds_id)
    idx <- idx1 & (splicings$exon_strand == "+")
    utr_end[idx] <- splicings$cds_start[idx] - 1L
    idx <- idx1 & (splicings$exon_strand == "-")
    utr_start[idx] <- splicings$cds_end[idx] + 1L

    ## split by grouping variable
    ans <- .makeUTRsByTranscript(txdb, splicings, utr_start, utr_end)
    ans <- .set_group_names(ans, use.names, txdb, "tx")
    ans <- .assignMetadataList(ans, txdb)
    ans
}

.make3UTRsByTranscript <- function(txdb, splicings, use.names=FALSE)
{
    exons_with_cds <- which(!is.na(splicings$cds_id))
    idx <- .exons_with_3utr(splicings$tx_id, exons_with_cds)
    splicings <- S4Vectors:::extract_data_frame_rows(splicings, idx)

    ## Compute the UTR starts/ends.
    utr_start <- splicings$exon_start
    utr_end <- splicings$exon_end
    idx1 <- !is.na(splicings$cds_id)
    idx <- idx1 & (splicings$exon_strand == "+")
    utr_start[idx] <- splicings$cds_end[idx] + 1L
    idx <- idx1 & (splicings$exon_strand == "-")
    utr_end[idx] <- splicings$cds_start[idx] - 1L

    ## split by grouping variable
    ans <- .makeUTRsByTranscript(txdb, splicings, utr_start, utr_end)
    ans <- .set_group_names(ans, use.names, txdb, "tx")
    ans <- .assignMetadataList(ans, txdb)
    ans
}

setGeneric("fiveUTRsByTranscript",
    function(x, ...) standardGeneric("fiveUTRsByTranscript")
)

setMethod("fiveUTRsByTranscript", "TxDb",
    function(x, use.names=FALSE)
    {
        splicings <- .getSplicingsForTranscriptsWithCDSs(x)
        .make5UTRsByTranscript(x, splicings, use.names=use.names)
    }
)

setGeneric("threeUTRsByTranscript",
    function(x, ...) standardGeneric("threeUTRsByTranscript")
)

setMethod("threeUTRsByTranscript", "TxDb",
    function(x, use.names=FALSE)
    {
        splicings <- .getSplicingsForTranscriptsWithCDSs(x)
        .make3UTRsByTranscript(x, splicings, use.names=use.names)
    }
)

