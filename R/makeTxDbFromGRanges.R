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
###   - required: type, ID
###   - optional: Parent, Name, Dbxref, geneID
### Used in R/makeTxDbFromGFF.R
GFF3_COLNAMES <- c("type", "ID", "Parent", "Name", "Dbxref", "geneID")

### Expected metadata columns for GRanges in GTF format:
###   - required: type, gene_id, transcript_id
###   - optional: exon_id
### Used in R/makeTxDbFromGFF.R
GTF_COLNAMES <- c("type", "gene_id", "transcript_id", "exon_id")

.GENE_TYPES <- c("gene", "pseudogene", "transposable_element_gene")
.TX_TYPES <- c("transcript", "pseudogenic_transcript", "primary_transcript",
               "mRNA", "ncRNA", "rRNA", "snoRNA", "snRNA", "tRNA", "tmRNA",
               "miRNA", "miRNA_primary_transcript",
               "RNase_P_RNA", "RNase_MRP_RNA", "SRP_RNA", "misc_RNA",
               "antisense_RNA", "antisense")
.EXON_TYPES <- c("exon", "pseudogenic_exon")
.CDS_TYPES <- "CDS"
.STOP_CODON_TYPES <- "stop_codon"

### Used in R/makeTxDbFromGFF.R
GFF_FEATURE_TYPES <- c(.GENE_TYPES, .TX_TYPES, .EXON_TYPES,
                       .CDS_TYPES, .STOP_CODON_TYPES)

.get_type <- function(gr_mcols)
{
    type <- gr_mcols$type
    if (is.null(type))
        stop("'gr' must have a \"type\" metadata column")
    ## Return a factor where all levels are in use.
    if (is.factor(type)) {
        levels_in_use <- levels(type)[tabulate(type, nbins=nlevels(type)) != 0L]
    } else {
        if (!is.character(type))
            stop(wmsg("the \"type\" metadata column must be ",
                      "a character vector or factor"))
        levels_in_use <- unique(type)
    }
    factor(type, levels=levels_in_use)
}

### Return a character vector or NULL.
.get_gene_id <- function(gr_mcols)
{
    gene_id <- gr_mcols$gene_id
    if (is.null(gene_id))
        return(NULL)
    if (!is.character(gene_id))
        gene_id <- as.character(gene_id)
    gene_id
}

.no_id <- function(id) {is.null(id) || all(is.na(id))}

### Return a character vector, or integer vector (if inferred ids), or NULL.
.get_transcript_id <- function(gr_mcols, gene_id, type)
{
    transcript_id <- gr_mcols$transcript_id
    ## We've seen silly GTF files that contain only lines of type transcript
    ## but no transcript_id tag.
    if (!.no_id(gene_id) && .no_id(transcript_id) && all(type %in% .TX_TYPES))
        return(seq_along(type))  # inferred ids
    if (is.null(transcript_id))
        return(NULL)
    if (!is.character(transcript_id))
        transcript_id <- as.character(transcript_id)
    transcript_id
}

### If we have no "ID" metadata column but "gene_id" and "transcript_id"
### metadata columns (and if they don't contain only NAs) then we assume the
### GRanges object is in GTF format. Otherwise we assume it's in GFF3 format.
.is_gtf_format <- function(ID, gene_id, transcript_id)
    (.no_id(ID) && !.no_id(gene_id) && !.no_id(transcript_id))

.get_ID <- function(gr_mcols, type, gene_id, transcript_id, gtf.format=FALSE)
{
    if (gtf.format) {
        ## GTF format
        ID <- rep.int(NA_character_, nrow(gr_mcols))
        idx1 <- which(type %in% .GENE_TYPES)
        ID[idx1] <- gene_id[idx1]
        idx2 <- which(type %in% .TX_TYPES)
        ID[idx2] <- transcript_id[idx2]
        exon_id <- gr_mcols$exon_id
        if (!is.null(exon_id)) {
            idx3 <- which(type %in% .EXON_TYPES)
            ID[idx3] <- exon_id[idx3]
        }
    } else {
        ## GFF3 format
        ID <- gr_mcols$ID
        if (is.null(ID))
            ID <- rep.int(NA_character_, nrow(gr_mcols))
        if (!is.character(ID))
            ID <- as.character(ID)
        ID[ID %in% ""] <- NA_character_
    }
    ID
}

.get_Parent <- function(gr_mcols, type, gene_id, transcript_id,
                        gtf.format=FALSE)
{
    if (gtf.format) {
        ## GTF format
        Parent <- character(nrow(gr_mcols))
        idx1 <- which(type %in% .TX_TYPES)
        Parent[idx1] <- gene_id[idx1]
        idx2 <- which(type %in% c(.EXON_TYPES, .CDS_TYPES, .STOP_CODON_TYPES))
        Parent[idx2] <- transcript_id[idx2]
    } else {
        ## GFF3 format
        Parent <- gr_mcols$Parent
        if (is.null(Parent))
            Parent <- rep.int(NA_character_, nrow(gr_mcols))
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

### Remove "<type_level>:" prefix (like in "gene:Si016158m.g" and
### "CDS:Si016158m") if present in *all* non-NA elements of 'ID'.
.infer_Name_from_ID <- function(ID, type_level)
{
    prefix <- paste0(type_level, ":")
    if (!all(substr(ID, start=1L, stop=nchar(prefix)) == prefix, na.rm=TRUE))
        return(ID)
    substr(ID, start=nchar(prefix)+1L, stop=nchar(ID))
}

.get_Name <- function(gr_mcols, type, ID, transcript_id)
{
    Name <- gr_mcols$Name
    if (is.null(Name))
        Name <- rep.int(NA_character_, nrow(gr_mcols))
    if (!is.character(Name))
        Name <- as.character(Name)
    Name[Name %in% ""] <- NA_character_
    ## If, for a given type (e.g. "exon"), none of the features of that type
    ## has a Name (i.e. Name is NA for all the features of that type), then we
    ## infer their Name from their ID.
    for (type_level in levels(type)) {
        ## We make an exception for transcripts with inferred ids.
        if ((type_level %in% .TX_TYPES) && is.integer(transcript_id))
            next
        idx <- which(type == type_level)
        if (all(is.na(Name[idx])))
            Name[idx] <- .infer_Name_from_ID(ID[idx], type_level)
    }
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

### The FlyBase people use the geneID tag to assign an external gene id to
### each of their transcripts in their GFF3 files. There is no corresponding
### "gene" line in the file.
.get_geneID <- function(gr_mcols)
{
    geneID <- gr_mcols$geneID
    if (is.null(geneID))
        return(NULL)
    if (!is.character(geneID))
        geneID <- as.character(geneID)
    geneID
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Get the gene, cds, stop_codon, exon, and transcript indices.
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

.get_stop_codon_IDX <- function(type)
{
    which(type %in% .STOP_CODON_TYPES)
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

### noextx :== transcript with no exons
.get_noextx_IDX <- function(tx_IDX, ID, exon_IDX, Parent)
{
    tx_IDX[!(ID[tx_IDX] %in% unlist(Parent[exon_IDX], use.names=FALSE))]
}

### CDS in transcripts with no exons will also be considered exons (i.e. added
### to the set of exons).
.get_cds_with_noextx_parent_IDX <- function(cds_IDX, Parent, noextx_IDX, ID)
{
    ## as.logical() will fail if a CDS has at least 2 parents and 1 is a
    ## transcript with no exons and the other is not.
    cds_IDX[as.logical(unique(Parent[cds_IDX] %in% ID[noextx_IDX]))]
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
### Extract the 'exons', 'cds', and 'stop_codons' data frames
###

.prepare_dropped_msg <- function(dropped, what)
{
    msg <- c("The following ", what, " were dropped")
    if (nrow(dropped) > 6L)
        msg <- c(msg, " (showing only the 6 first)")
    c(msg, ":\n", paste(capture.output(print(head(dropped, n=6L))),
                                       collapse="\n"))
}

.extract_rank_from_id <- function(id, parent_id)
{
    if (length(id) == 0L)
        return(integer(0))
    id_parts <- strsplit(id, "\\.|:")
    ## Fix non-sensical output of strsplit() on empty strings.
    id_parts[elementNROWS(id_parts) == 0L] <- ""
    unlisted_id_parts <- unlist(id_parts, use.names=FALSE)
    idx <- cumsum(elementNROWS(id_parts))
    rank <- unlisted_id_parts[idx]
    rank <- suppressWarnings(as.integer(rank))
    if (any(is.na(rank)))
        return(NULL)
    ## Sanity check.
    tmp <- unname(splitAsList(rank, parent_id))
    if (any(any(duplicated(tmp))))
        return(NULL)
    if (!(all(min(tmp) == 1L) && all(max(tmp) == elementNROWS(tmp))))
        return(NULL)
    rank
}

.infer_rank_from_position <- function(tx_id, exon_chrom, exon_strand,
                                             exon_start, exon_end)
{
    chrom_by_tx <- split(Rle(exon_chrom), tx_id)
    tx_chrom <- runValue(chrom_by_tx)
    bad_tx1 <- names(which(elementNROWS(tx_chrom) > 1L))

    strand_by_tx <- split(Rle(exon_strand), tx_id)
    tx_strand <- runValue(strand_by_tx)
    is_bad <- elementNROWS(tx_strand) > 1L
    bad_tx2 <- names(which(is_bad))
    tx_strand[is_bad] <- "*"
    minus_idx <- which(as.character(tx_strand) == "-")

    ex_by_tx <- split(IRanges(exon_start, exon_end), tx_id)
    reduced_ex_by_tx <- reduce(ex_by_tx, min.gapwidth=0L)
    bad_tx3 <- names(which(elementNROWS(reduced_ex_by_tx) !=
                           elementNROWS(ex_by_tx)))

    start_by_tx <- start(ex_by_tx)
    start_by_tx[minus_idx] <- start_by_tx[minus_idx] * (-1L)
    rank <- .rank.CompressedList(start_by_tx, ties.method="first")
    ans <- unsplit(rank, tx_id)

    bad_tx <- unique(c(bad_tx1, bad_tx2, bad_tx3))
    ans[tx_id %in% bad_tx] <- NA_integer_

    ans
}

### Can be used to extract exons, cds, or stop codons.
.extract_exons_from_GRanges <- function(exon_IDX, gr, ID, Name, Parent,
                                        feature=c("exon", "cds", "stop_codon"),
                                        gtf.format=FALSE)
{
    feature <- match.arg(feature)
    what <- switch(feature, exon="exon", cds="CDS", stop_codon="stop codon")
    exon_Parent <- Parent[exon_IDX]
    if (any(any(duplicated(exon_Parent))))
        stop(wmsg("some ", what, "s are mapped twice to the same transcript"))

    tx_id <- factor(unlist(exon_Parent, use.names=FALSE))
    nparent_per_ex <- elementNROWS(exon_Parent)
    exon_id <- rep.int(ID[exon_IDX], nparent_per_ex)
    exon_name <- rep.int(Name[exon_IDX], nparent_per_ex)
    exon_chrom <- rep.int(seqnames(gr)[exon_IDX], nparent_per_ex)
    exon_strand <- rep.int(strand(gr)[exon_IDX], nparent_per_ex)
    exon_start <- rep.int(start(gr)[exon_IDX], nparent_per_ex)
    exon_end <- rep.int(end(gr)[exon_IDX], nparent_per_ex)

    if (!gtf.format) {
        ## Drop orphan exons (or orphan cds or stop codons).
        is_orphan <- !(tx_id %in% ID)
        norphan <- sum(is_orphan)
        if (norphan != 0L) {
            dropped <- data.frame(seqid=exon_chrom[is_orphan],
                                  start=exon_start[is_orphan],
                                  end=exon_end[is_orphan],
                                  strand=exon_strand[is_orphan],
                                  ID=exon_id[is_orphan],
                                  Parent=tx_id[is_orphan],
                                  Name=exon_name[is_orphan],
                                  stringsAsFactors=FALSE)
            warning(.prepare_dropped_msg(dropped, paste("orphan", what)))
            keep_idx <- which(!is_orphan)
            tx_id <- tx_id[keep_idx]
            exon_id <- exon_id[keep_idx]
            exon_name <- exon_name[keep_idx]
            exon_chrom <- exon_chrom[keep_idx]
            exon_strand <- exon_strand[keep_idx]
            exon_start <- exon_start[keep_idx]
            exon_end <- exon_end[keep_idx]
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

    if (feature == "exon") {
        exon_rank <- .extract_rank_from_id(exon_id, tx_id)
        if (is.null(exon_rank))
            exon_rank <- .infer_rank_from_position(tx_id,
                                                   exon_chrom, exon_strand,
                                                   exon_start, exon_end)
        exons$exon_rank <- exon_rank
    } else {
        colnames(exons) <- sub("^exon_", "", colnames(exons))
    }
    exons
}

### For each transcript with no exons, infer a single exon that is the same as
### the transcript itself. Other heuristics might have been used earlier to
### give exons to these poor transcripts with no exons (like inferring the
### exons from the CDS, to support the GeneDB case) and .add_missing_exons()
### should only be used as the last and desperate solution to give them an
### exon after the other heuristics have failed to give them one.
.add_missing_exons <- function(exons, transcripts)
{
    is_noextx <- !(transcripts$tx_id %in% exons$tx_id)
    missing_exons <- S4Vectors:::extract_data_frame_rows(transcripts, is_noextx)
    if (nrow(missing_exons) == 0L)
        return(exons)
    missing_exons$tx_type <- NULL
    TX2EXON_COLNAMES <- c(
        tx_id="tx_id",
        tx_name="exon_name",
        tx_chrom="exon_chrom",
        tx_strand="exon_strand",
        tx_start="exon_start",
        tx_end="exon_end"
    )
    colnames(missing_exons) <- TX2EXON_COLNAMES[colnames(missing_exons)]
    missing_exons$exon_rank <- 1L
    rbind(exons, missing_exons)
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
    bad_tx <- names(which(elementNROWS(tx_name) != 1L))
    if (length(bad_tx) != 0L) {
        in1string <- paste0(sort(bad_tx), collapse=", ")
        stop(wmsg("The following transcripts have multiple parts that cannot ",
                  "be merged because of incompatible Name: ", in1string))
    }
    tx_name <- as.character(tx_name)

    tx_type <- unique(transcripts_by_id[ , "tx_type"])
    bad_tx <- names(which(elementNROWS(tx_type) != 1L))
    if (length(bad_tx) != 0L) {
        in1string <- paste0(sort(bad_tx), collapse=", ")
        stop(wmsg("The following transcripts have multiple parts that cannot ",
                  "be merged because of incompatible type: ", in1string))
    }
    tx_type <- as.character(tx_type)

    tx_chrom <- unique(transcripts_by_id[ , "tx_chrom"])
    bad_tx <- names(which(elementNROWS(tx_chrom) != 1L))
    if (length(bad_tx) != 0L) {
        in1string <- paste0(sort(bad_tx), collapse=", ")
        stop(wmsg("The following transcripts have multiple parts that cannot ",
                  "be merged because of incompatible seqnames: ", in1string))
    }
    tx_chrom <- as.character(tx_chrom)

    tx_strand <- unique(transcripts_by_id[ , "tx_strand"])
    bad_tx <- names(which(elementNROWS(tx_strand) != 1L))
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
        tx_type=tx_type,
        tx_chrom=tx_chrom,
        tx_strand=tx_strand,
        tx_start=tx_start,
        tx_end=tx_end,
        stringsAsFactors=FALSE
    )
}

.extract_transcripts_from_GRanges <- function(tx_IDX, gr, type, ID, Name)
{
    tx_id <- ID[tx_IDX]
    transcripts <- data.frame(
        tx_id=tx_id,
        tx_name=Name[tx_IDX],
        tx_type=type[tx_IDX],
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
    ## We've seen silly GFF/GTF files where each transcript in the file is
    ## reported to be on its own contig and spans it (start=1) but no strand
    ## is reported for the transcript. We set the strand to "+" for all these
    ## transcripts.
    if (all(transcripts$tx_strand == "*")
     && anyDuplicated(transcripts$tx_chrom) == 0L
     && all(transcripts$tx_start == 1L)) {
        transcripts$tx_strand <- rep.int(strand("+"), nrow(transcripts))
    }
    transcripts
}

.infer_transcripts_from_exons <- function(exons)
{
    exons_tx_id <- exons$tx_id
    if (!is.character(exons_tx_id))
        exons_tx_id <- as.character(exons_tx_id)
    if (!is.character(exons$exon_chrom))
        exons$exon_chrom <- as.character(exons$exon_chrom)
    if (!is.character(exons$exon_strand))
        exons$exon_strand <- as.character(exons$exon_strand)
    exons_by_id <- splitAsList(exons, exons_tx_id, drop=TRUE)

    tx_id <- names(exons_by_id)

    tx_type <- rep.int("inferred_from_exons", length(tx_id))

    tx_chrom <- unique(exons_by_id[ , "exon_chrom"])
    tx_chrom[elementNROWS(tx_chrom) != 1L] <- NA_character_
    tx_chrom <- as.character(tx_chrom)

    tx_strand <- unique(exons_by_id[ , "exon_strand"])
    tx_strand[elementNROWS(tx_strand) != 1L] <- NA_character_
    tx_strand <- as.character(tx_strand)

    tx_start <- unname(min(exons_by_id[ , "exon_start"]))
    tx_end <- unname(max(exons_by_id[ , "exon_end"]))
    
    data.frame(
        tx_id=tx_id,
        tx_name=tx_id,
        tx_type=tx_type,
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
    orphan_exons <- S4Vectors:::extract_data_frame_rows(exons, is_orphan)
    missing_transcripts <- .infer_transcripts_from_exons(orphan_exons)
    rbind(transcripts, missing_transcripts)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Make the 'splicings' data frame
###

.find_exon_cds <- function(exons, cds, what="CDS")
{
    query <- GRanges(cds$chrom,
                     IRanges(cds$start, cds$end),
                     cds$strand)
    subject <- GRanges(exons$exon_chrom,
                       IRanges(exons$exon_start, exons$exon_end),
                       exons$exon_strand)
    hits <- findOverlaps(query, subject, type="within")
    hits <- hits[cds$tx_id[queryHits(hits)] == exons$tx_id[subjectHits(hits)]]
    q_hits <- queryHits(hits)
    s_hits <- subjectHits(hits)

    bad_tx <- unique(cds$tx_id[q_hits[duplicated(q_hits)]])
    if (length(bad_tx) != 0L) {
        because <- c("they have ", what, "s that are mapped to ",
                     "more than one exon")
        .reject_transcripts(bad_tx, because)
    }

    cds2exon <- selectHits(hits, select="arbitrary")
    bad_tx <- unique(cds$tx_id[is.na(cds2exon)])
    if (length(bad_tx) != 0L) {
        because <- c("they have ", what, "s that cannot be mapped to an exon")
        .reject_transcripts(bad_tx, because)
    }

    bad_tx <- unique(exons$tx_id[s_hits[duplicated(s_hits)]])
    if (length(bad_tx) != 0L) {
        in1string <- paste0(sort(bad_tx), collapse=", ")
        warning(wmsg("The following transcripts have exons that contain ",
                     "more than one ", what, " (only the first ", what,
                     " was kept for each exon): ", in1string))
    }
    selectHits(t(hits), select="first")
}

.make_splicings <- function(exons, cds, stop_codons=NULL)
{
    cds_tx_id <- factor(cds$tx_id, levels=levels(exons$tx_id))
    if (any(is.na(cds_tx_id)))
        stop(wmsg("some CDS cannot be mapped to an exon"))
    cds$tx_id <- cds_tx_id

    exon2cds <- .find_exon_cds(exons, cds)
    cds_name <- cds$name[exon2cds]
    cds_start <- cds$start[exon2cds]
    cds_end <- cds$end[exon2cds]

    if (!is.null(stop_codons)) {
        stop_codons_tx_id <- factor(stop_codons$tx_id,
                                    levels=levels(exons$tx_id))
        if (any(is.na(stop_codons_tx_id)))
            stop(wmsg("some stop codons cannot be mapped to an exon"))
        stop_codons$tx_id <- stop_codons_tx_id

        exon2stop_codon <- .find_exon_cds(exons, stop_codons,
                                          what="stop codon")
        stop_codon_name <- stop_codons$name[exon2stop_codon]
        stop_codon_start <- stop_codons$start[exon2stop_codon]
        stop_codon_end <- stop_codons$end[exon2stop_codon]

        ## Exons with no CDS get the stop codon as CDS.
        replace_idx <- which(is.na(exon2cds) & !is.na(exon2stop_codon))
        cds_name[replace_idx] <- stop_codon_name[replace_idx]
        cds_start[replace_idx] <- stop_codon_start[replace_idx]
        cds_end[replace_idx] <- stop_codon_end[replace_idx]

        ## Exons with a CDS and a stop codon have the latter merged into the
        ## former.
        merge_idx <- which(!is.na(exon2cds) & !is.na(exon2stop_codon))
        start1 <- cds_start[merge_idx]
        end1 <- cds_end[merge_idx]
        start2 <- stop_codon_start[merge_idx]
        end2 <- stop_codon_end[merge_idx]
        gap <- pmax.int(start1, start2) - pmin.int(end1, end2) - 1L
        bad_tx <- unique(exons$tx_id[merge_idx[gap > 0L]])
        if (length(bad_tx) != 0L) {
            because <- "they have incompatible CDS and stop codons"
            .reject_transcripts(bad_tx, because)
        }
        cds_start[merge_idx] <- pmin.int(start1, start2)
        cds_end[merge_idx] <- pmax.int(end1, end2)
    }
    cbind(exons, cds_name, cds_start, cds_end, stringsAsFactors=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the 'genes' data frame
###

.extract_genes_from_gff3_GRanges <- function(gene_IDX, tx_IDX,
                                             ID, Name, Parent,
                                             Dbxref=NULL, geneID=NULL)
{
    tx_id <- ID[tx_IDX]
    tx2genes <- Parent[tx_IDX]

    ## Genes that we consider transcripts are therefore their own parents.
    idx0 <- which(tx_IDX %in% gene_IDX)
    tx2genes[idx0] <- tx_id[idx0]

    ## Transcripts with no parents are sometimes linked to a gene via the
    ## Dbxref or geneID tag.

    ## First we try to find their parent via the Dbxref tag, if present.
    if (!is.null(Dbxref)) {
        idx0 <- which(elementNROWS(tx2genes) == 0L)
        tx_Dbxref <- Dbxref[tx_IDX[idx0]]
        gene_Dbxref <- Dbxref[gene_IDX]
        tx_Dbxref_unlisted <- unlist(tx_Dbxref, use.names=FALSE)
        gene_Dbxref_unlisted <- unlist(gene_Dbxref, use.names=FALSE)
        hits <- findMatches(tx_Dbxref_unlisted, gene_Dbxref_unlisted)
        hits <- remapHits(hits, Lnodes.remapping=togroup(tx_Dbxref),
                                new.nLnode=length(tx_Dbxref),
                                Rnodes.remapping=togroup(gene_Dbxref),
                                new.nRnode=length(gene_Dbxref))
        tx2genes[idx0] <- relist(ID[gene_IDX][subjectHits(hits)],
                                 as(hits, "PartitioningByEnd"))
    }

    ## Replace gene IDs by gene Names and remove NAs.
    tmp <- unlist(tx2genes, use.names=FALSE)
    ID2Name <- Name[gene_IDX]
    names(ID2Name) <- ID[gene_IDX]
    tx2genes <- relist(unname(ID2Name[tmp]), tx2genes)
    tx2genes <- tx2genes[!is.na(tx2genes)]

    ## Then, if we still have transcripts with no parent, we use the geneID
    ## tag (if present) to assign them an external gene id.
    if (!is.null(geneID)) {
        idx0 <- which(elementNROWS(tx2genes) == 0L)
        tx2genes[idx0] <- geneID[tx_IDX[idx0]]
        tx2genes <- tx2genes[!is.na(tx2genes)]
    }

    tx_id <- rep.int(tx_id, elementNROWS(tx2genes))
    gene_id <- unlist(tx2genes, use.names=FALSE)
    data.frame(tx_id=tx_id, gene_id=gene_id, stringsAsFactors=FALSE)
}

.extract_genes_from_gtf_GRanges <- function(transcript_id, gene_id,
                                            transcripts)
{
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

.normarg_metadata <- function(metadata, Genome, taxonomyId)
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
    if (!is.na(taxonomyId)) {
        GenomeInfoDb:::.checkForAValidTaxonomyId(taxonomyId)
        org <- GenomeInfoDb:::.lookupSpeciesFromTaxId(taxonomyId)
        organism <- paste(org[1,2], org[1,3])
        idx <- match("Organism", metadata$name)
        if (!is.na(idx))
            metadata <- S4Vectors:::extract_data_frame_rows(metadata, -idx)
        idx <- match("Taxonomy ID", metadata$name)
        if (!is.na(idx))
            metadata <- S4Vectors:::extract_data_frame_rows(metadata, -idx)
        df2 <- data.frame(name=c("Organism", "Taxonomy ID"),
                          value=c(organism, taxonomyId),
                          stringsAsFactors=FALSE)
        metadata <- rbind(metadata, df2)
    }
    metadata
}

### If 'drop.stop.codons' is TRUE then the "stop_codon" lines are ignored.
### Otherwise (the default) the stop codons are considered to be part of the
### CDS and merged to them.
makeTxDbFromGRanges <- function(gr, drop.stop.codons=FALSE, metadata=NULL,
                                taxonomyId=NA)
{
    if (!is(gr, "GenomicRanges"))
        stop("'gr' must be a GRanges object")

    Genome <- unique(genome(gr))
    if (length(Genome) > 1L)
        stop("all the sequences in 'seqinfo(gr)' must belong ",
             "to the same genome")
    if (length(Genome) == 0L)
        Genome <- NA_character_

    if (!isTRUEorFALSE(drop.stop.codons))
        stop("'drop.stop.codons' must be TRUE or FALSE")

    metadata <- .normarg_metadata(metadata, Genome, taxonomyId)

    gr_mcols <- mcols(gr)

    ## Get the metadata columns of interest.
    type <- .get_type(gr_mcols)
    gene_id <- .get_gene_id(gr_mcols)
    transcript_id <- .get_transcript_id(gr_mcols, gene_id, type)
    gtf.format <- .is_gtf_format(gr_mcols$ID, gene_id, transcript_id)
    ID <- .get_ID(gr_mcols, type, gene_id, transcript_id,
                  gtf.format=gtf.format)
    Parent <- .get_Parent(gr_mcols, type, gene_id, transcript_id,
                          gtf.format=gtf.format)
    Name <- .get_Name(gr_mcols, type, ID, transcript_id)

    ## Get the gene, cds, stop_codon, exon, and transcript indices.
    gene_IDX <- .get_gene_IDX(type)
    cds_IDX <- .get_cds_IDX(type)
    cds_with_gene_parent_IDX <- .get_cds_with_gene_parent_IDX(cds_IDX,
                                          Parent, gene_IDX, ID,
                                          gtf.format=gtf.format)
    if (!drop.stop.codons)
        stop_codon_IDX <- .get_stop_codon_IDX(type)
    exon_IDX <- .get_exon_IDX(type, cds_with_gene_parent_IDX)
    exon_with_gene_parent_IDX <- .get_exon_with_gene_parent_IDX(exon_IDX,
                                          Parent, gene_IDX, ID,
                                          gtf.format=gtf.format)
    gene_as_tx_IDX <- .get_gene_as_tx_IDX(gene_IDX, ID,
                                          exon_with_gene_parent_IDX, Parent)
    tx_IDX <- .get_tx_IDX(type, gene_as_tx_IDX)

    noextx_IDX <- .get_noextx_IDX(tx_IDX, ID, exon_IDX, Parent)
    if (length(noextx_IDX) != 0L) {
        ## Infer exons for transcripts with no exons (noextx).

        ## For a noextx with CDS (the GeneDB case): infer 1 exon per CDS that
        ## is the same as the CDS itself.
        cds_with_noextx_parent_IDX <- .get_cds_with_noextx_parent_IDX(
                                              cds_IDX, Parent,
                                              noextx_IDX, ID)
        exon_IDX <- sort(c(exon_IDX, cds_with_noextx_parent_IDX))

        ## For a noextx with no CDS (the miRBase case): we'll take care of
        ## them later with .add_missing_exons().
    }

    ## Extract the 'exons', 'cds', 'stop_codons', 'transcripts',
    ## and 'genes' data frames.
    exons <- .extract_exons_from_GRanges(exon_IDX, gr, ID, Name, Parent,
                                   feature="exon", gtf.format=gtf.format)
    cds <- .extract_exons_from_GRanges(cds_IDX, gr, ID, Name, Parent,
                                   feature="cds", gtf.format=gtf.format)
    if (!drop.stop.codons) {
        stop_codons <- .extract_exons_from_GRanges(stop_codon_IDX,
                                   gr, ID, Name, Parent,
                                   feature="stop_codon", gtf.format=gtf.format)
    } else {
        stop_codons <- NULL
    }
    transcripts <- .extract_transcripts_from_GRanges(tx_IDX, gr,
                                                     type, ID, Name)
    if (gtf.format) {
        transcripts <- .add_missing_transcripts(transcripts, exons)
        genes <- .extract_genes_from_gtf_GRanges(transcript_id, gene_id,
                                                 transcripts)
    } else {
        Dbxref <- .get_Dbxref(gr_mcols)
        geneID <- .get_geneID(gr_mcols)
        genes <- .extract_genes_from_gff3_GRanges(gene_IDX, tx_IDX,
                                                  ID, Name, Parent,
                                                  Dbxref, geneID)
    }

    ## Drop hopeless transcripts.
    drop1_tx <- as.character(exons$tx_id[is.na(exons$exon_rank)])
    drop2_tx <- as.character(transcripts$tx_id[is.na(transcripts$tx_chrom)])
    drop3_tx <- as.character(transcripts$tx_id[is.na(transcripts$tx_strand)])
    drop1_tx <- unique(drop1_tx)
    drop2_tx <- unique(drop2_tx)
    drop3_tx <- unique(drop3_tx)
    drop_tx <- unique(c(drop1_tx, drop2_tx, drop3_tx))
    if (length(drop_tx) != 0L) {
        exons <- S4Vectors:::extract_data_frame_rows(exons,
                                     !(exons$tx_id %in% drop_tx))
        cds <- S4Vectors:::extract_data_frame_rows(cds,
                                     !(cds$tx_id %in% drop_tx))
        stop_codons <- S4Vectors:::extract_data_frame_rows(stop_codons,
                                     !(stop_codons$tx_id %in% drop_tx))
        transcripts <- S4Vectors:::extract_data_frame_rows(transcripts,
                                     !(transcripts$tx_id %in% drop_tx))
        genes <- S4Vectors:::extract_data_frame_rows(genes,
                                     !(genes$tx_id %in% drop_tx))
        if (length(drop1_tx) != 0L) {
            in1string <- paste0(sort(drop1_tx), collapse=", ")
            warning(wmsg(
                "The following transcripts were dropped because their exon ",
                "ranks could not be inferred (either because the exons are ",
                "not on the same chromosome/strand or because they are not ",
                "separated by introns): ", in1string))
        }
        if (length(drop2_tx) != 0L) {
            in1string <- paste0(sort(drop2_tx), collapse=", ")
            warning(wmsg(
                "The following transcripts were dropped because no genomic ",
                "ranges could be found for them and their ranges could not ",
                "be inferred from their exons either (because they have them ",
                "on more than one chromosome): ", in1string))
        }
        if (length(drop3_tx) != 0L) {
            in1string <- paste0(sort(drop3_tx), collapse=", ")
            warning(wmsg(
                "The following transcripts were dropped because no genomic ",
                "ranges could be found for them and their ranges could not ",
                "be inferred from their exons either (because they have them ",
                "on both strands): ", in1string))
        }
    }

    ## We might still have transcripts with no exons (noextx). Take care of
    ## them now.
    exons <- .add_missing_exons(exons, transcripts)

    ## .make_splicings() can also reject some transcripts.
    .flush_rejected_tx_envir()
    splicings <- .make_splicings(exons, cds, stop_codons)
    drop_tx <- .get_rejected_transcripts()
    if (length(drop_tx) != 0L) {
        transcripts <- S4Vectors:::extract_data_frame_rows(transcripts,
                                           !(transcripts$tx_id %in% drop_tx))
        splicings <- S4Vectors:::extract_data_frame_rows(splicings,
                                         !(splicings$tx_id %in% drop_tx))
        genes <- S4Vectors:::extract_data_frame_rows(genes,
                                     !(genes$tx_id %in% drop_tx))
    }

    ## Turn the "tx_id" column in 'splicings', 'genes', and 'transcripts' into
    ## an integer vector.
    ## TODO: Maybe makeTxDb() could take care of that.
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some non-exported tools to test makeTxDbFromGRanges() on real data.
###
### NOTE: These tests take too long to be part of the formal unit tests!
###

.light_txdb_dump <- function(txdb)
{
    txdb_dump <- as.list(txdb)

    keep_cols <- c("tx_name", "tx_chrom", "tx_strand", "tx_start", "tx_end")
    transcripts <- txdb_dump$transcripts[keep_cols]
    oo <- order(transcripts$tx_chrom,
                transcripts$tx_start, transcripts$tx_end,
                transcripts$tx_strand,
                transcripts$tx_name)
    transcripts <- S4Vectors:::extract_data_frame_rows(transcripts, oo)

    keep_cols <- c("exon_rank", "exon_name", "exon_chrom", "exon_strand",
                   "exon_start", "exon_end", "cds_start", "cds_end")
    splicings <- txdb_dump$splicings[keep_cols]
    oo <- order(splicings$exon_chrom,
                splicings$exon_start, splicings$exon_end,
                splicings$exon_strand,
                splicings$exon_name)
    splicings <- S4Vectors:::extract_data_frame_rows(splicings, oo)

    list(transcripts=transcripts, splicings=splicings)
}

test_makeTxDbFromGRanges_on_Ensembl_organism_gtf <- function(organism)
{
    cat("Download GTF file ... ")
    gtf_url <- ftp_url_to_Ensembl_gtf()
    url <- paste0(gtf_url, organism, "/")
    organism_files <- ls_ftp_url(url)
    gtf_file <- grep("gtf", organism_files, ignore.case=TRUE, value=TRUE)
    if (length(gtf_file) != 1L) {
        cat("ABORT (0 or more than 1 GTF file found)\n")
        return(NA)
    }
    url <- paste0(url, gtf_file)
    local_gtf_file <- file.path(tempdir(), gtf_file)
    download.file(url, local_gtf_file)

    cat("Import ", local_gtf_file, " as GRanges object 'gr' ... ", sep="")
    gr <- import(local_gtf_file, colnames=GTF_COLNAMES,
                                 feature.type=GFF_FEATURE_TYPES)
    cat("\n")

    cat("txdb1 <- makeTxDbFromGRanges(gr) ... ", sep="")
    txdb1 <- makeTxDbFromGRanges(gr)
    cat("\n")

    dataset <- strsplit(organism, "_", fixed=TRUE)[[1L]]
    dataset <- paste0(substr(dataset[[1L]], 1L, 1L),
                      dataset[[2L]],
                      "_gene_ensembl")
    cat("txdb2 <- makeTxDbFromBiomart(\"", dataset, "\") ... ", sep="")
    txdb2 <- try(makeTxDbFromBiomart(dataset))
    if (is(txdb2, "try-error")) {
        cat("ABORT\n")
        return(NA)
    }
    cat("\n")

    cat("'txdb1' and 'txdb2' have the same transcripts, exons, and CDS: ")
    dump1 <- .light_txdb_dump(txdb1)
    dump2 <- .light_txdb_dump(txdb2)
    OK <- identical(dump1, dump2)
    if (OK) {
        cat("YES\n")
    } else {
        cat("NO!!!!!!\n")
    }
    return(OK)
}

test_makeTxDbFromGRanges_on_Ensembl_gtf <- function(all=FALSE)
{
    gtf_url <- ftp_url_to_Ensembl_gtf()
    if (all) {
        organisms <- ls_ftp_url(gtf_url)
    } else {
        organisms <- c("caenorhabditis_elegans",
                       "ciona_intestinalis",
                       "danio_rerio",
                       "drosophila_melanogaster",
                       "gallus_gallus",
                       "homo_sapiens",
                       "mus_musculus",
                       "rattus_norvegicus",
                       "saccharomyces_cerevisiae")
    }
    norganism <- length(organisms)
    for (i in seq_len(norganism)) {
        organism <- organisms[[i]]
        cat("\n")
        cat("*************************************************************\n")
        cat("[", i , "/", norganism, "] ",
            "Testing makeTxDbFromGRanges() on ", organism, "\n", sep="")
        cat("*************************************************************\n")
        cat("\n")
        test_makeTxDbFromGRanges_on_Ensembl_organism_gtf(organism)
    }
    cat("\n")
    cat("DONE.\n")
}


if (FALSE) {
library(GenomicFeatures)
source("GenomicFeatures/R/makeTxDbFromGRanges.R")
library(rtracklayer)

## Test with GRanges obtained from GFF3 files
## ==========================================

GFF3_files <- system.file("extdata", "GFF3_files", package="GenomicFeatures")

file1 <- file.path(GFF3_files, "TheCanonicalGene_v1.gff3")
gr1 <- import(file1, format="gff3", colnames=GFF3_COLNAMES,
                     feature.type=GFF_FEATURE_TYPES)
txdb1 <- makeTxDbFromGRanges(gr1)
txdb1

file2 <- file.path(GFF3_files, "TheCanonicalGene_v2.gff3")
gr2 <- import(file2, format="gff3", colnames=GFF3_COLNAMES,
                     feature.type=GFF_FEATURE_TYPES)
txdb2 <- makeTxDbFromGRanges(gr2)
txdb2

file3 <- file.path(GFF3_files, "a.gff3")
gr3 <- import(file3, format="gff3", colnames=GFF3_COLNAMES,
                     feature.type=GFF_FEATURE_TYPES)
txdb3 <- makeTxDbFromGRanges(gr3)
txdb3

file4 <- file.path(GFF3_files, "dmel-1000-r5.11.filtered.gff")
gr4 <- import(file4, format="gff3", colnames=GFF3_COLNAMES,
                     feature.type=GFF_FEATURE_TYPES)
txdb4 <- makeTxDbFromGRanges(gr4)
txdb4  # exactly the same as with makeTxDbFromGFF()

file5 <- file.path(GFF3_files, "NC_011025.gff")
gr5 <- import(file5, format="gff3", colnames=GFF3_COLNAMES,
                     feature.type=GFF_FEATURE_TYPES)
txdb5 <- makeTxDbFromGRanges(gr5)
txdb5

genome <- "GCA_000364345.1_Macaca_fascicularis_5.0"
filename <- paste0(genome, "_genomic.gff.gz")
url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/", genome, "/", filename)
file6 <- file.path(tempdir(), file6)
download.file(url, file6)
gr6 <- import(file6, format="gff3", colnames=GFF3_COLNAMES,
                     feature.type=GFF_FEATURE_TYPES)

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
gr1 <- import(file1, format="gtf", colnames=GTF_COLNAMES,
                     feature.type=GFF_FEATURE_TYPES)
txdb1 <- makeTxDbFromGRanges(gr1)
txdb1

file2 <- file.path(GTF_files, "Aedes_aegypti.partial.gtf")
gr2 <- import(file2, format="gtf", colnames=GTF_COLNAMES,
                     feature.type=GFF_FEATURE_TYPES)
txdb2 <- makeTxDbFromGRanges(gr2)
txdb2
}

