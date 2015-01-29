### =========================================================================
### makeTxDbFromGRanges()
### -------------------------------------------------------------------------
###


### TODO: Move this to IRanges and make it the S4 "rank" method for
### CompressedList objects. Note that the same trick can be used to implement
### a fast "rank" method for list and SimpleList objects. So it would be
### great to implement a .rank.List() in S4Vectors that does the same as
### below but without using IRanges idioms like PartitioningByEnd() and
### togroup().
.rank.CompressedList <- function(x, na.last=TRUE,
                                 ties.method=c("average", "first",
                                               "random", "max", "min"))
{
    ties.method <- match.arg(ties.method)
    if (!identical(ties.method, "first"))
        stop(wmsg("\"rank\" method for CompressedList objects ",
                  "only supports 'ties.method=\"first\"' at the moment"))
    x_partitioning <- PartitioningByEnd(x)
    unlisted_x <- unlist(x, use.names=FALSE)
    oo <- S4Vectors:::orderIntegerPairs(togroup(x_partitioning), unlisted_x)
    unlisted_ans <- integer(length(oo))
    unlisted_ans[oo] <- seq_along(oo)
    ans <- relist(unlisted_ans, x_partitioning)
    x_len <- length(x)
    if (x_len != 0L) {
        offsets <- c(0L, end(x_partitioning)[-x_len])
        ans <- ans - offsets
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get the metadata columns of interest
###
### Expected metadata columns for GRanges in GFF3 format:
###   - required: type, ID, Parent
###   - optional: Name, Dbxref
###
### Expected metadata columns for GRanges in GTF format:
###   - type, gene_id, transcript_id, exon_id
###

.GENE_TYPES <- c("gene", "pseudogene")
.TX_TYPES <- c("mRNA", "transcript", "primary_transcript",
               "ncRNA", "rRNA", "snoRNA", "snRNA", "tRNA", "tmRNA")
.EXON_TYPES <- "exon"
.CDS_TYPES <- "CDS"

.get_type <- function(gr_mcols)
{
    type <- gr_mcols$type
    if (is.null(type))
        stop("'gr' must have a \"type\" metadata column")
    if (!is.factor(type))
        type <- factor(type, levels=unique(type))
    type
}

.get_gene_id <- function(gr_mcols)
{
    gene_id <- gr_mcols$gene_id
    if (is.null(gene_id))
        return(NULL)
    if (!is.character(gene_id))
        gene_id <- as.character(gene_id)
    gene_id
}

.get_transcript_id <- function(gr_mcols)
{
    transcript_id <- gr_mcols$transcript_id
    if (is.null(transcript_id))
        return(NULL)
    if (!is.character(transcript_id))
        transcript_id <- as.character(transcript_id)
    transcript_id
}

.get_ID <- function(gr_mcols, type, gene_id, transcript_id)
{
    ID <- gr_mcols$ID
    if (is.null(ID)) {
        ## If there's no ID column then we assume the GRanges object is in
        ## GTF format.
        if (is.null(gene_id) || is.null(transcript_id))
            stop("'gr' must have an \"ID\" metadata column, ",
                 "or a \"gene_id\" and \"transcript_id\" metadata column")
        ID <- character(nrow(gr_mcols))
        idx1 <- which(type %in% .GENE_TYPES)
        ID[idx1] <- gene_id[idx1]
        idx2 <- which(type %in% .TX_TYPES)
        ID[idx2] <- transcript_id[idx2]
        exon_id <- gr_mcols$exon_id
        if (!is.null(exon_id)) {
            idx3 <- which(type %in% .EXON_TYPES)
            ID[idx3] <- exon_id[idx3]
        }
        return(ID)
    }
    if (!is.character(ID))
        ID <- as.character(ID)
    ID
}

### If there's no Parent column then we assume the GRanges object is in
### GTF format.
.is_gtf_format <- function(gr_mcols)
    is.null(gr_mcols$Parent)

.get_Parent <- function(gr_mcols, type, gene_id, transcript_id,
                        gtf.format=FALSE)
{
    if (gtf.format) {
        if (is.null(gene_id) || is.null(transcript_id))
            stop("'gr' must have a \"Parent\" metadata column, ",
                 "or a \"gene_id\" and \"transcript_id\" metadata column")
        Parent <- character(nrow(gr_mcols))
        idx1 <- which(type %in% .TX_TYPES)
        Parent[idx1] <- gene_id[idx1]
        idx2 <- which(type %in% c(.EXON_TYPES, .CDS_TYPES))
        Parent[idx2] <- transcript_id[idx2]
    } else {
        Parent <- gr_mcols$Parent
    }
    if (!is(Parent, "CompressedCharacterList")) {
        Parent0 <- Parent
        Parent <- as(Parent, "CompressedCharacterList")
        if (is.atomic(Parent0)) {
            na_idx <- which(is.na(Parent0))
            Parent[na_idx] <- CharacterList(character(0))
        }
    }
    Parent
}

.get_Name <- function(gr_mcols, ID)
{
    Name <- gr_mcols$Name
    if (is.null(Name))
        Name <- ID
    if (!is.character(Name))
        Name <- as.character(Name)
    Name[Name %in% ""] <- NA_character_
    Name
}

.get_Dbxref <- function(gr_mcols)
{
    Dbxref <- gr_mcols$Dbxref
    if (is.null(Dbxref))
        return(NULL)
    if (is(Dbxref, "List") || is.list(Dbxref)) {
        if (!is(Dbxref, "CompressedCharacterList"))
            Dbxref <- as(Dbxref, "CompressedCharacterList")
    } else {
        if (!is.character(Dbxref))
            Dbxref <- as.character(Dbxref)
    }
    Dbxref
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get the transcript, exon, cds, and gene indices
###

.get_gene_IDX <- function(type)
{
    which(type %in% .GENE_TYPES)
}

.get_cds_IDX <- function(type)
{
    which(type %in% .CDS_TYPES)
}

### Returns the index of CDS whose Parent is a gene (this happens with some
### GFF files e.g. inst/extdata/GFF3_files/NC_011025.gff). These CDS will
### also be considered exons (i.e. added to the set of exons).
.get_cds_with_gene_parent_IDX <- function(cds_IDX, Parent, gene_IDX, ID,
                                          gtf.format=FALSE)
{
    if (gtf.format)
        return(integer(0))
    cds_IDX[as.logical(unique(Parent[cds_IDX] %in% ID[gene_IDX]))]
}

.get_exon_IDX <- function(type, cds_with_gene_parent_IDX)
{
    sort(c(which(type %in% .EXON_TYPES), cds_with_gene_parent_IDX))
}

### Returns the index of exons whose Parent is a gene.
.get_exon_with_gene_parent_IDX <- function(exon_IDX, Parent, gene_IDX, ID,
                                           gtf.format=FALSE)
{
    if (gtf.format)
        return(integer(0))
    exon_IDX[as.logical(unique(Parent[exon_IDX] %in% ID[gene_IDX]))]
}

### Returns the index of genes that have exons as direct children. These genes
### will also be considered transcripts (i.e. added to the set of transcripts).
.get_gene_as_tx_IDX <- function(gene_IDX, ID, exon_with_gene_parent_IDX, Parent)
{
    gene_IDX[ID[gene_IDX] %in% unlist(Parent[exon_with_gene_parent_IDX],
                                      use.names=FALSE)]
}

.get_tx_IDX <- function(type, gene_as_tx_IDX)
{
    sort(c(which(type %in% .TX_TYPES), gene_as_tx_IDX))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Management of rejected transcripts
###

.rejected_tx_envir <- new.env(hash=TRUE, parent=emptyenv())

.ls_rejected_tx_envir <- function()
    ls(envir=.rejected_tx_envir, all.names=TRUE, sorted=FALSE)

.flush_rejected_tx_envir <- function()
{
    objnames <- .ls_rejected_tx_envir()
    rm(list=objnames, envir=.rejected_tx_envir)
}

.reject_transcripts <- function(tx_ids, because)
{
    tx_ids <- sort(as.character(tx_ids))
    in1string <- paste0(tx_ids, collapse=", ")
    warning(wmsg("The following transcripts were rejected because ",
                 because, ": ", in1string))

    objnames <- .ls_rejected_tx_envir()
    nobj <- length(objnames)
    if (nobj == 0L) {
        new_objname <- 1L
    } else {
        new_objname <- as.integer(objnames[nobj]) + 1L
    }
    new_objname <- sprintf("%08d", new_objname)
    assign(new_objname, tx_ids, envir=.rejected_tx_envir)
}

.get_rejected_transcripts <- function()
{
    objnames <- .ls_rejected_tx_envir()
    unique(unlist(mget(objnames, envir=.rejected_tx_envir), use.names=FALSE))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the 'exons' and 'cds' data frames
###

.extract_rank_from_id <- function(id, parent_id)
{
    if (length(id) == 0L)
        return(integer(0))
    id_parts <- strsplit(id, "\\.|:")
    ## Fix non-sensical output of strsplit() on empty strings.
    id_parts[elementLengths(id_parts) == 0L] <- ""
    unlisted_id_parts <- unlist(id_parts, use.names=FALSE)
    idx <- cumsum(elementLengths(id_parts))
    rank <- unlisted_id_parts[idx]
    rank <- suppressWarnings(as.integer(rank))
    if (any(is.na(rank)))
        return(NULL)
    ## Sanity check.
    tmp <- unname(splitAsList(rank, parent_id))
    if (any(any(duplicated(tmp))))
        return(NULL)
    if (!(all(min(tmp) == 1L) && all(max(tmp) == elementLengths(tmp))))
        return(NULL)
    rank
}

.infer_rank_from_position <- function(tx_id, exon_chrom, exon_strand,
                                             exon_start, exon_end)
{
    chrom_by_tx <- split(Rle(exon_chrom), tx_id)
    tx_chrom <- runValue(chrom_by_tx)
    bad_tx1 <- names(which(elementLengths(tx_chrom) > 1L))

    strand_by_tx <- split(Rle(exon_strand), tx_id)
    tx_strand <- runValue(strand_by_tx)
    is_bad <- elementLengths(tx_strand) > 1L
    bad_tx2 <- names(which(is_bad))
    tx_strand[is_bad] <- "*"
    minus_idx <- which(as.character(tx_strand) == "-")

    ex_by_tx <- split(IRanges(exon_start, exon_end), tx_id)
    reduced_ex_by_tx <- reduce(ex_by_tx, min.gapwidth=0L)
    bad_tx3 <- names(which(elementLengths(reduced_ex_by_tx) !=
                           elementLengths(ex_by_tx)))

    start_by_tx <- start(ex_by_tx)
    start_by_tx[minus_idx] <- start_by_tx[minus_idx] * (-1L)
    rank <- .rank.CompressedList(start_by_tx, ties.method="first")
    ans <- unsplit(rank, tx_id)

    bad_tx <- unique(c(bad_tx1, bad_tx2, bad_tx3))
    ans[tx_id %in% bad_tx] <- NA_integer_

    ans
}

### Can be used to extract exons or cds.
.extract_exons_from_GRanges <- function(exon_IDX, gr, ID, Name, Parent,
                                        gtf.format=FALSE, for.cds=FALSE)
{
    feature.type <- if (for.cds) "CDS" else "exon"
    exon_Parent <- Parent[exon_IDX]
    if (any(any(duplicated(exon_Parent))))
        stop(wmsg("some ", feature.type, "s are mapped twice ",
                  "to the same transcript"))

    tx_id <- factor(unlist(exon_Parent, use.names=FALSE))
    nparent_per_ex <- elementLengths(exon_Parent)
    exon_id <- rep.int(ID[exon_IDX], nparent_per_ex)
    exon_name <- rep.int(Name[exon_IDX], nparent_per_ex)
    exon_chrom <- rep.int(seqnames(gr)[exon_IDX], nparent_per_ex)
    exon_strand <- rep.int(strand(gr)[exon_IDX], nparent_per_ex)
    exon_start <- rep.int(start(gr)[exon_IDX], nparent_per_ex)
    exon_end <- rep.int(end(gr)[exon_IDX], nparent_per_ex)

    if (!gtf.format) {
        ## Drop orphan exons (or cds).
        is_orphan <- !(tx_id %in% ID)
        norphan <- sum(is_orphan)
        if (norphan != 0L) {
            keep_idx <- which(!is_orphan)
            tx_id <- tx_id[keep_idx]
            exon_id <- exon_id[keep_idx]
            exon_name <- exon_name[keep_idx]
            exon_chrom <- exon_chrom[keep_idx]
            exon_strand <- exon_strand[keep_idx]
            exon_start <- exon_start[keep_idx]
            exon_end <- exon_end[keep_idx]
            warning(wmsg(norphan, " orphan ", feature.type, "s were dropped"))
        }
    }

    exons <- data.frame(
        tx_id=tx_id,
        exon_name=exon_name,
        exon_chrom=exon_chrom,
        exon_strand=exon_strand,
        exon_start=exon_start,
        exon_end=exon_end,
        stringsAsFactors=FALSE
    )

    if (for.cds) {
        colnames(exons) <- sub("^exon", "cds", colnames(exons))
    } else {
        exon_rank <- .extract_rank_from_id(exon_id, tx_id)
        if (is.null(exon_rank))
            exon_rank <- .infer_rank_from_position(tx_id,
                                                   exon_chrom, exon_strand,
                                                   exon_start, exon_end)
        exons$exon_rank <- exon_rank
    }
    exons
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the 'transcripts' data frame
###

.merge_transcript_parts <- function(transcripts)
{
    tx_id <- transcripts$tx_id
    if (!is.character(tx_id))
        tx_id <- as.character(tx_id)
    transcripts_by_id <- splitAsList(transcripts, tx_id, drop=TRUE)

    tx_id <- names(transcripts_by_id)

    tx_name <- unique(transcripts_by_id[ , "tx_name"])
    bad_tx <- names(which(elementLengths(tx_name) != 1L))
    if (length(bad_tx) != 0L) {
        in1string <- paste0(sort(bad_tx), collapse=", ")
        stop(wmsg("The following transcripts have multiple parts that cannot ",
                  "be merged because of incompatible Name: ", in1string))
    }
    tx_name <- as.character(tx_name)

    tx_chrom <- unique(transcripts_by_id[ , "tx_chrom"])
    bad_tx <- names(which(elementLengths(tx_chrom) != 1L))
    if (length(bad_tx) != 0L) {
        in1string <- paste0(sort(bad_tx), collapse=", ")
        stop(wmsg("The following transcripts have multiple parts that cannot ",
                  "be merged because of incompatible seqnames: ", in1string))
    }
    tx_chrom <- as.character(tx_chrom)

    tx_strand <- unique(transcripts_by_id[ , "tx_strand"])
    bad_tx <- names(which(elementLengths(tx_strand) != 1L))
    if (length(bad_tx) != 0L) {
        in1string <- paste0(sort(bad_tx), collapse=", ")
        stop(wmsg("The following transcripts have multiple parts that cannot ",
                  "be merged because of incompatible strand: ", in1string))
    }
    tx_strand <- as.character(tx_strand)

    tx_start <- unname(min(transcripts_by_id[ , "tx_start"]))
    tx_end <- unname(max(transcripts_by_id[ , "tx_end"]))
    
    data.frame(
        tx_id=tx_id,
        tx_name=tx_name,
        tx_chrom=tx_chrom,
        tx_strand=tx_strand,
        tx_start=tx_start,
        tx_end=tx_end,
        stringsAsFactors=FALSE
    )
}

.extract_transcripts_from_GRanges <- function(tx_IDX, gr, ID, Name)
{
    tx_id <- ID[tx_IDX]
    transcripts <- data.frame(
        tx_id=tx_id,
        tx_name=Name[tx_IDX],
        tx_chrom=seqnames(gr)[tx_IDX],
        tx_strand=strand(gr)[tx_IDX],
        tx_start=start(gr)[tx_IDX],
        tx_end=end(gr)[tx_IDX],
        stringsAsFactors=FALSE
    )
    merged_tx <- unique(tx_id[duplicated(tx_id)])
    if (length(merged_tx) != 0L) {
        transcripts <- .merge_transcript_parts(transcripts)
        in1string <- paste0(sort(merged_tx), collapse=", ")
        warning(wmsg("The following transcripts have multiple parts that ",
                     "were merged: ", in1string))
    }
    transcripts
}

.infer_transcripts_from_exons <- function(exons)
{
    exons_tx_id <- exons$tx_id
    if (!is.character(exons_tx_id))
        exons_tx_id <- as.character(exons_tx_id)
    exons_by_id <- splitAsList(exons, exons_tx_id, drop=TRUE)

    tx_id <- names(exons_by_id)

    tx_chrom <- unique(exons_by_id[ , "exon_chrom"])
    bad_tx <- names(which(elementLengths(tx_chrom) != 1L))
    if (length(bad_tx) != 0L) {
        in1string <- paste0(sort(bad_tx), collapse=", ")
        stop(wmsg("No genomic ranges found for the following transcripts ",
                  "and the ranges could not be inferred from their exons ",
                  "(because the transcripts have exons on more than one ",
                  "chromosome): ", in1string))
    }
    tx_chrom <- as.character(tx_chrom)

    tx_strand <- unique(exons_by_id[ , "exon_strand"])
    bad_tx <- names(which(elementLengths(tx_strand) != 1L))
    if (length(bad_tx) != 0L) {
        in1string <- paste0(sort(bad_tx), collapse=", ")
        stop(wmsg("No genomic ranges found for the following transcripts ",
                  "and the ranges could not be inferred from their exons ",
                  "(because the transcripts have exons on both strands): ",
                  in1string))
    }
    tx_strand <- as.character(tx_strand)

    tx_start <- unname(min(exons_by_id[ , "exon_start"]))
    tx_end <- unname(max(exons_by_id[ , "exon_end"]))
    
    data.frame(
        tx_id=tx_id,
        tx_name=tx_id,
        tx_chrom=tx_chrom,
        tx_strand=tx_strand,
        tx_start=tx_start,
        tx_end=tx_end,
        stringsAsFactors=FALSE
    )
}

.add_missing_transcripts <- function(transcripts, exons)
{
    is_orphan <- !(exons$tx_id %in% transcripts$tx_id)
    orphan_exons <- exons[is_orphan, ]
    missing_transcripts <- .infer_transcripts_from_exons(orphan_exons)
    rbind(transcripts, missing_transcripts)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Make the 'splicings' data frame
###

.find_exon_cds <- function(exons, cds)
{
    query <- GRanges(cds$cds_chrom,
                     IRanges(cds$cds_start, cds$cds_end),
                     cds$cds_strand)
    subject <- GRanges(exons$exon_chrom,
                       IRanges(exons$exon_start, exons$exon_end),
                       exons$exon_strand)
    hits <- findOverlaps(query, subject, type="within")
    hits <- hits[cds$tx_id[queryHits(hits)] == exons$tx_id[subjectHits(hits)]]
    q_hits <- queryHits(hits)
    s_hits <- subjectHits(hits)

    bad_tx <- unique(cds$tx_id[q_hits[duplicated(q_hits)]])
    if (length(bad_tx) != 0L) {
        because <- "they have CDS that are mapped to more than one exon"
        .reject_transcripts(bad_tx, because)
    }

    cds2exon <- selectHits(hits, select="arbitrary")
    bad_tx <- unique(cds$tx_id[is.na(cds2exon)])
    if (length(bad_tx) != 0L) {
        because <- "they have CDS that cannot be mapped to an exon"
        .reject_transcripts(bad_tx, because)
    }

    bad_tx <- unique(exons$tx_id[s_hits[duplicated(s_hits)]])
    if (length(bad_tx) != 0L) {
        in1string <- paste0(sort(bad_tx), collapse=", ")
        warning(wmsg("The following transcripts have exons that contain ",
                     "more than one CDS (only the first CDS was kept ",
                     "for each exon): ", in1string))
    }
    selectHits(t(hits), select="first")
}

.make_splicings <- function(exons, cds)
{
    cds_tx_id <- factor(cds$tx_id, levels=levels(exons$tx_id))
    if (any(is.na(cds_tx_id)))
        stop(wmsg("some CDS cannot be mapped to an exon"))
    cds$tx_id <- cds_tx_id

    exon2cds <- .find_exon_cds(exons, cds)
    cds_name <- cds$cds_name[exon2cds]
    cds_start <- cds$cds_start[exon2cds]
    cds_end <- cds$cds_end[exon2cds]
    cbind(exons, cds_name, cds_start, cds_end, stringsAsFactors=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the 'genes' data frame
###

.extract_genes_from_gff3_GRanges <- function(gene_IDX, tx_IDX, ID, Name,
                                             Parent, Dbxref)
{
    tx_id <- ID[tx_IDX]
    tx2genes <- Parent[tx_IDX]

    ## Genes that we consider transcripts are therefore their own parents.
    idx0 <- which(tx_IDX %in% gene_IDX)
    tx2genes[idx0] <- tx_id[idx0]

    ## Transcripts with no parents are sometimes linked to a gene via the
    ## Dbxref tag.
    if (!is.null(Dbxref)) {
        idx0 <- which(elementLengths(tx2genes) == 0L)
        tx_Dbxref <- Dbxref[tx_IDX][idx0]
        gene_Dbxref <- Dbxref[gene_IDX]
        tx_Dbxref_unlisted <- unlist(tx_Dbxref, use.names=FALSE)
        gene_Dbxref_unlisted <- unlist(gene_Dbxref, use.names=FALSE)
        hits <- findMatches(tx_Dbxref_unlisted, gene_Dbxref_unlisted)
        hits <- remapHits(hits, query.map=togroup(tx_Dbxref),
                                new.queryLength=length(tx_Dbxref),
                                subject.map=togroup(gene_Dbxref),
                                new.subjectLength=length(gene_Dbxref))
        tx2genes[idx0] <- relist(ID[gene_IDX][subjectHits(hits)],
                                 as(hits, "PartitioningByEnd"))
    }

    tx_id <- rep.int(tx_id, elementLengths(tx2genes))
    geneID2Name <- Name[gene_IDX]
    names(geneID2Name) <- ID[gene_IDX]
    gene_id <- unname(geneID2Name[unlist(tx2genes, use.names=FALSE)])
    keep_idx <- which(!is.na(gene_id))
    tx_id <- tx_id[keep_idx]
    gene_id <- gene_id[keep_idx]
    data.frame(tx_id=tx_id, gene_id=gene_id, stringsAsFactors=FALSE)
}

.extract_genes_from_gtf_GRanges <- function(transcript_id, gene_id,
                                            transcripts)
{
    if (is.null(transcript_id) || is.null(gene_id)) {
        transcript_id <- character(0)
        gene_id <- character(0)
    }

    ## Keep only unique (tx_id, gene_id) pairs.
    tx_id <- factor(transcript_id, levels=unique(transcript_id))
    gene_id <- factor(gene_id, levels=unique(gene_id))
    keep_idx <- which(!S4Vectors:::duplicatedIntegerPairs(tx_id, gene_id))
    tx_id <- tx_id[keep_idx]
    gene_id <- gene_id[keep_idx]

    ## Keep only (tx_id, gene_id) pairs that point to a row in 'transcripts'.
    keep_idx <- which(tx_id %in% transcripts$tx_id)
    tx_id <- tx_id[keep_idx]
    gene_id <- gene_id[keep_idx]

    data.frame(tx_id=tx_id, gene_id=gene_id, stringsAsFactors=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the 'chrominfo' data frame
###

.extract_chrominfo_from_GRanges <- function(gr)
{
    gr_seqinfo <- seqinfo(gr)
    data.frame(
        chrom=seqnames(gr_seqinfo),
        length=seqlengths(gr_seqinfo),
        is_circular=isCircular(gr_seqinfo),
        stringsAsFactors=FALSE
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTxDbFromGRanges()
###

.normarg_metadata <- function(metadata, Genome)
{
    if (!is.null(metadata)) {
        if (!is.data.frame(metadata))
            stop("'metadata' must be NULL or a data.frame")
        if (!setequal(colnames(metadata), c("name", "value")))
            stop("'metadata' columns must be \"name\" and \"value\"")
    }
    if (!("Genome" %in% metadata$name)) {
        df1 <- data.frame(name="Genome", value=Genome, stringsAsFactors=FALSE)
        metadata <- rbind(metadata, df1)
    }
    metadata
}

makeTxDbFromGRanges <- function(gr, metadata=NULL)
{
    if (!is(gr, "GenomicRanges"))
        stop("'gr' must be a GRanges object")

    Genome <- unique(genome(gr))
    if (length(Genome) != 1L)
        stop("all the sequences in 'seqinfo(gr)' must belong ",
             "to the same genome")
    metadata <- .normarg_metadata(metadata, Genome)

    gr_mcols <- mcols(gr)

    ## Get the metadata columns of interest.
    type <- .get_type(gr_mcols)
    gene_id <- .get_gene_id(gr_mcols)
    transcript_id <- .get_transcript_id(gr_mcols)
    ID <- .get_ID(gr_mcols, type, gene_id, transcript_id)
    gtf.format <- .is_gtf_format(gr_mcols)
    Parent <- .get_Parent(gr_mcols, type, gene_id, transcript_id,
                          gtf.format=gtf.format)
    Name <- .get_Name(gr_mcols, ID)
    Dbxref <- .get_Dbxref(gr_mcols)

    ## Get the transcript, exon, cds, and gene indices.
    gene_IDX <- .get_gene_IDX(type)
    cds_IDX <- .get_cds_IDX(type)
    cds_with_gene_parent_IDX <- .get_cds_with_gene_parent_IDX(cds_IDX,
                                          Parent, gene_IDX, ID,
                                          gtf.format=gtf.format)
    exon_IDX <- .get_exon_IDX(type, cds_with_gene_parent_IDX)
    exon_with_gene_parent_IDX <- .get_exon_with_gene_parent_IDX(exon_IDX,
                                          Parent, gene_IDX, ID,
                                          gtf.format=gtf.format)
    gene_as_tx_IDX <- .get_gene_as_tx_IDX(gene_IDX, ID,
                                          exon_with_gene_parent_IDX, Parent)
    tx_IDX <- .get_tx_IDX(type, gene_as_tx_IDX)

    ## Extract the 'exons', 'cds', 'transcripts', and 'genes' data frames.
    exons <- .extract_exons_from_GRanges(exon_IDX, gr, ID, Name, Parent,
                                         gtf.format=gtf.format)
    cds <- .extract_exons_from_GRanges(cds_IDX, gr, ID, Name, Parent,
                                         gtf.format=gtf.format, for.cds=TRUE)
    transcripts <- .extract_transcripts_from_GRanges(tx_IDX, gr, ID, Name)
    if (gtf.format) {
        transcripts <- .add_missing_transcripts(transcripts, exons)
        genes <- .extract_genes_from_gtf_GRanges(transcript_id, gene_id,
                                                 transcripts)
    } else {
        genes <- .extract_genes_from_gff3_GRanges(gene_IDX, tx_IDX,
                                                  ID, Name, Parent, Dbxref)
    }

    ## Drop transcripts for which the exon ranks could not be inferred.
    drop_tx <- unique(as.character(exons$tx_id[is.na(exons$exon_rank)]))
    if (length(drop_tx) != 0L) {
        exons <- exons[!(exons$tx_id %in% drop_tx), ]
        cds <- cds[!(cds$tx_id %in% drop_tx), ]
        transcripts <- transcripts[!(transcripts$tx_id %in% drop_tx), ]
        genes <- genes[!(genes$tx_id %in% drop_tx), ]
        in1string <- paste0(sort(drop_tx), collapse=", ")
        warning(wmsg(
            "The following transcripts were dropped because their exon ",
            "ranks could not be inferred (either because the exons are ",
            "not on the same chromosome/strand or because they are not ",
            "separated by introns): ", in1string))
    }

    .flush_rejected_tx_envir()
    splicings <- .make_splicings(exons, cds)
    drop_tx <- .get_rejected_transcripts()
    if (length(drop_tx) != 0L) {
        transcripts <- transcripts[!(transcripts$tx_id %in% drop_tx), ]
        splicings <- splicings[!(splicings$tx_id %in% drop_tx), ]
        genes <- genes[!(genes$tx_id %in% drop_tx), ]
    }

    ## Turn the "tx_id" column in 'splicings', 'genes', and 'transcripts' into
    ## an integer vector. TODO: Maybe makeTxDb() could take care of this.
    splicings_tx_id <- match(splicings$tx_id, transcripts$tx_id)
    if (any(is.na(splicings_tx_id)))
        stop(wmsg("some exons are linked to transcripts ",
                  "not found in the file"))
    splicings$tx_id <- splicings_tx_id
    genes_tx_id <- match(genes$tx_id, transcripts$tx_id)
    if (any(is.na(genes_tx_id)))
        stop(wmsg("some genes are linked to transcripts ",
                  "not found in the file"))
    genes$tx_id <- genes_tx_id
    transcripts$tx_id <- seq_along(transcripts$tx_id)

    chrominfo <- .extract_chrominfo_from_GRanges(gr)

    makeTxDb(transcripts, splicings,
             genes=genes, chrominfo=chrominfo,
             metadata=metadata, reassign.ids=TRUE)
}


if (FALSE) {
library(GenomicFeatures)
source("GenomicFeatures/R/makeTxDbFromGRanges.R")
feature.type <- c(.GENE_TYPES, .TX_TYPES, .EXON_TYPES, .CDS_TYPES)
library(rtracklayer)

## Test with GRanges obtained from GFF3 files
## ==========================================

GFF3_files <- system.file("extdata", "GFF3_files", package="GenomicFeatures")

file1 <- file.path(GFF3_files, "TheCanonicalGene_v1.gff3")
gr1 <- import(file1, format="gff3", feature.type=feature.type)
txdb1 <- makeTxDbFromGRanges(gr1)
txdb1

file2 <- file.path(GFF3_files, "TheCanonicalGene_v2.gff3")
gr2 <- import(file2, format="gff3", feature.type=feature.type)
txdb2 <- makeTxDbFromGRanges(gr2)
txdb2

file3 <- file.path(GFF3_files, "a.gff3")
gr3 <- import(file3, format="gff3", feature.type=feature.type)
txdb3 <- makeTxDbFromGRanges(gr3)
txdb3

file4 <- file.path(GFF3_files, "dmel-1000-r5.11.filtered.gff")
gr4 <- import(file4, format="gff3", feature.type=feature.type)
txdb4 <- makeTxDbFromGRanges(gr4)
txdb4  # exactly the same as with makeTxDbFromGFF()

file5 <- file.path(GFF3_files, "NC_011025.gff")
gr5 <- import(file5, format="gff3", feature.type=feature.type)
txdb5 <- makeTxDbFromGRanges(gr5)
txdb5

genome <- "GCA_000364345.1_Macaca_fascicularis_5.0"
filename <- paste0(genome, "_genomic.gff.gz")
url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/", genome, "/", filename)
file6 <- file.path(tempdir(), file6)
download.file(url, file6)
gr6 <- import(file6, format="gff3", feature.type=feature.type)

## Compared with makeTxDbFromGFF():
##
## (a) makeTxDbFromGFF() fails on TheCanonicalGene_v1.gff3
##
## (b) gene_id, tx_name, exon_name, and cds_name are now imported from the
##     Name tag instead of the ID tag (GFF3 Spec: "IDs do not have meaning
##     outside the file in which they reside")


## Test with GRanges obtained from GTF files
## =========================================

GTF_files <- system.file("extdata", "GTF_files", package="GenomicFeatures")

## test1.gtf grabbed from http://mblab.wustl.edu/GTF22.html (5 exon gene with
## 3 translated exons).
file1 <- file.path(GTF_files, "test1.gtf")
gr1 <- import(file1, format="gtf", feature.type=feature.type)
txdb1 <- makeTxDbFromGRanges(gr1)
txdb1

file2 <- file.path(GTF_files, "Aedes_aegypti.partial.gtf")
gr2 <- import(file2, format="gtf", feature.type=feature.type)
txdb2 <- makeTxDbFromGRanges(gr2)
txdb2
}

