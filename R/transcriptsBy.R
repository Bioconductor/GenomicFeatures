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
    grl <- keep_user_seqlevels_from_TxDb(grl, txdb)
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
    S4Vectors:::fancy_mseq(lengths, offset=offset)
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
    offset <- last_exon_with_cds - 1L
    lengths <- cumsum(nexon) - offset
    S4Vectors:::fancy_mseq(lengths, offset=offset)
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
    keep_user_seqlevels_from_TxDb(grl, txdb)
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

