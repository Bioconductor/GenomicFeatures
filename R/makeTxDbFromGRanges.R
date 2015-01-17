### =========================================================================
### makeTxDbFromGRanges()
### -------------------------------------------------------------------------
###


.GENE_TYPES <- "gene"
.TX_TYPES <- c("mRNA", "ncRNA", "rRNA", "snoRNA", "snRNA", "tRNA", "tmRNA")
.EXON_TYPES <- "exon"
.CDS_TYPES <- "CDS"

.get_gene_IDX <- function(type)
{
    which(type %in% .GENE_TYPES)
}

.get_cds_IDX <- function(type)
{
    which(type %in% .CDS_TYPES)
}

### Returns the index of CDS that are direct children of genes (this happens
### with some GFF files e.g. inst/extdata/GFF3_files/NC_011025.gff). These CDS
### will also be considered exons (i.e. added to the set of exons).
.get_cds_as_exon_IDX <- function(cds_IDX, Parent, gene_IDX, ID)
{
    cds_IDX[as.logical(unique(Parent[cds_IDX] %in% ID[gene_IDX]))]
}

.get_exon_IDX <- function(type, cds_as_exon_IDX)
{
    sort(c(which(type %in% .EXON_TYPES), cds_as_exon_IDX))
}

### Returns the index of genes that have CDS as direct children. These genes
### will also be considered transcripts (i.e. added to the set of transcripts).
.get_gene_as_tx_IDX <- function(gene_IDX, ID, cds_as_exon_IDX, Parent)
{
    gene_IDX[ID[gene_IDX] %in% unlist(Parent[cds_as_exon_IDX], use.names=FALSE)]
}

.get_tx_IDX <- function(type, gene_as_tx_IDX)
{
    sort(c(which(type %in% .TX_TYPES), gene_as_tx_IDX))
}

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
### Extract the 'transcripts' data frame
###

.extract_transcripts_from_GRanges <- function(tx_IDX, gr, ID, Name)
{
    tx_id <- ID[tx_IDX]
    if (anyDuplicated(tx_id))
        stop(wmsg("transcript IDs are not unique"))
    data.frame(
        tx_id=tx_id,
        tx_name=Name[tx_IDX],
        tx_chrom=seqnames(gr)[tx_IDX],
        tx_strand=strand(gr)[tx_IDX],
        tx_start=start(gr)[tx_IDX],
        tx_end=end(gr)[tx_IDX],
        stringsAsFactors=FALSE
    )
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the 'splicings' data frame
###

.extract_rank_from_id <- function(id, parent_id, for.cds=FALSE)
{
    id_parts <- strsplit(id, "\\.|:")
    idx <- cumsum(elementLengths(id_parts))
    rank <- unlist(id_parts, use.names=FALSE)[idx]
    rank <- suppressWarnings(as.integer(rank))
    if (any(is.na(rank)))
        return(NULL)
    ## Sanity check.
    tmp <- unname(splitAsList(rank, parent_id))
    if (any(any(duplicated(tmp))))
        return(NULL)
    if (for.cds)
        return(rank)
    if (!(all(min(tmp) == 1L) && all(max(tmp) == elementLengths(tmp))))
        return(NULL)
    rank
}

.infer_rank_from_position <- function(tx_id, exon_chrom, exon_strand,
                                             exon_start, exon_end)
{
    chrom_per_tx <- split(Rle(exon_chrom), tx_id)
    tx_chrom <- runValue(chrom_per_tx)
    bad_tx1 <- unique(names(which(elementLengths(tx_chrom) > 1L)))

    strand_per_tx <- split(Rle(exon_strand), tx_id)
    tx_strand <- runValue(strand_per_tx)
    is_bad <- elementLengths(tx_strand) > 1L
    bad_tx2 <- unique(names(which(is_bad)))
    tx_strand[is_bad] <- "*"
    minus_idx <- which(as.character(tx_strand) == "-")

    ex_per_tx <- split(IRanges(exon_start, exon_end), tx_id)
    bad_tx3 <- unique(names(which(elementLengths(reduce(ex_per_tx)) !=
                                  elementLengths(ex_per_tx))))

    bad_tx <- unique(c(bad_tx1, bad_tx2, bad_tx3))
    if (length(bad_tx)) {
        bad_tx <- paste0(bad_tx, collapse=", ")
        warning(wmsg(
            "The following transcripts were dropped because their exon ",
            "ranks could not be inferred (either because the exons are ",
            "not on the same chromosome/strand or because they are not ",
            "separated by introns): ", sort(bad_tx)))
    }

    start_per_tx <- start(ex_per_tx)
    start_per_tx[minus_idx] <- start_per_tx[minus_idx] * (-1L)
    rank <- .rank.CompressedList(start_per_tx, ties.method="first")
    ans <- unsplit(rank, tx_id)
    ans[tx_id %in% bad_tx] <- NA_integer_
    ans
}

### Can be used to extract exons or cds.
.extract_exons <- function(exon_IDX, gr, ID, Name, Parent, for.cds=FALSE)
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

    ans <- data.frame(
        tx_id=tx_id,
        exon_name=exon_name,
        exon_chrom=exon_chrom,
        exon_strand=exon_strand,
        exon_start=exon_start,
        exon_end=exon_end,
        stringsAsFactors=FALSE
    )
    if (for.cds) {
        colnames(ans) <- sub("^exon", "cds", colnames(ans))
    } else {
        exon_rank <- .extract_rank_from_id(exon_id, tx_id, for.cds=for.cds)
        if (is.null(exon_rank))
            exon_rank <- .infer_rank_from_position(tx_id,
                                                   exon_chrom, exon_strand,
                                                   exon_start, exon_end)
        ans$exon_rank <- exon_rank
    }
    ans
}

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
    if (length(bad_tx)) {
        bad_tx <- paste0(bad_tx, collapse=", ")
        stop(wmsg("The following transcripts have CDS that are mapped ",
                  "to more than one exon: ", sort(bad_tx)))
    }

    cds2exon <- selectHits(hits, select="arbitrary")
    bad_tx <- unique(cds$tx_id[is.na(cds2exon)])
    if (length(bad_tx)) {
        bad_tx <- paste0(bad_tx, collapse=", ")
        stop(wmsg("The following transcripts have CDS that cannot ",
                  "be mapped to an exon: ", sort(bad_tx)))
    }

    bad_tx <- unique(exons$tx_id[s_hits[duplicated(s_hits)]])
    if (length(bad_tx)) {
        bad_tx <- paste0(bad_tx, collapse=", ")
        warning(wmsg("The following transcripts have exons containing ",
                     "more than one CDS (only the first CDS was kept ",
                     "for each exon): ", sort(bad_tx)))
    }
    selectHits(t(hits), select="first")
}

.extract_splicings_from_GRanges <- function(exon_IDX, cds_IDX,
                                            gr, ID, Name, Parent)
{
    exons <- .extract_exons(exon_IDX, gr, ID, Name, Parent)
    cds <- .extract_exons(cds_IDX, gr, ID, Name, Parent, for.cds=TRUE)

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

.extract_genes_from_GRanges <- function(gene_IDX, tx_IDX, ID, Name,
                                        Parent, Dbxref)
{
    geneID2Name <- Name[gene_IDX]
    names(geneID2Name) <- ID[gene_IDX]

    tx2genes <- Parent[tx_IDX]

    ## Genes that we consider transcripts are therefore their own parents.
    idx0 <- which(tx_IDX %in% gene_IDX)
    tx2genes[idx0] <- ID[tx_IDX][idx0]

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

    tx_id <- rep.int(factor(ID[tx_IDX]), elementLengths(tx2genes))
    gene_id <- unname(geneID2Name[unlist(tx2genes, use.names=FALSE)])
    data.frame(
        tx_id=tx_id,
        gene_id=gene_id,
        stringsAsFactors=FALSE
    )
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

### Tags we care about:
### o If GRanges obtained from GFF3: type, ID, Name, Parent, Dbxref
### o If GRanges obtained from GTF:
###     type, gene_id, transcript_id, exon_number, gene_name, transcript_name
makeTxDbFromGRanges <- function(gr, metadata=NULL)
{
    gr_mcols <- mcols(gr)
    type <- gr_mcols[ , "type"]
    ID <- gr_mcols[ , "ID"]
    Name <- gr_mcols$Name  # optional tag
    if (is.null(Name))
        Name <- ID
    Parent <- gr_mcols[ , "Parent"]
    Dbxref <- gr_mcols$Dbxref  # optional tag

    gene_IDX <- .get_gene_IDX(type)
    cds_IDX <- .get_cds_IDX(type)
    cds_as_exon_IDX <- .get_cds_as_exon_IDX(cds_IDX, Parent, gene_IDX, ID)
    exon_IDX <- .get_exon_IDX(type, cds_as_exon_IDX)
    gene_as_tx_IDX <- .get_gene_as_tx_IDX(gene_IDX, ID, cds_as_exon_IDX, Parent)
    tx_IDX <- .get_tx_IDX(type, gene_as_tx_IDX)

    transcripts <- .extract_transcripts_from_GRanges(tx_IDX, gr, ID, Name)
    splicings <- .extract_splicings_from_GRanges(exon_IDX, cds_IDX,
                                                 gr, ID, Name, Parent)
    genes <- .extract_genes_from_GRanges(gene_IDX, tx_IDX,
                                         ID, Name, Parent, Dbxref)
    chrominfo <- .extract_chrominfo_from_GRanges(gr)

    ## Turn 'tx_id' into an integer vector.
    ## TODO: makeTxDb() should take care of this.
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

    makeTxDb(transcripts, splicings,
             genes=genes, chrominfo=chrominfo,
             metadata=metadata, reassign.ids=TRUE)
}


if (FALSE) {
library(GenomicFeatures)
source("GenomicFeatures/R/makeTxDbFromGRanges.R")

library(rtracklayer)

## Test with GRanges obtained from GFF3 files
## ==========================================

feature.type <- c(.GENE_TYPES, .TX_TYPES, .EXON_TYPES, .CDS_TYPES)

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

## Compared with makeTxDbFromGFF():
##
## (a) makeTxDbFromGFF() fails on TheCanonicalGene_v1.gff3
##
## (b) gene_id, tx_name, exon_name, and cds_name are now imported from the
##     Name tag instead of the ID tag (GFF3 Spec: "IDs do not have meaning
##     outside the file in which they reside")


## Test with GRanges obtained from GTF files
## =========================================

file <- system.file("extdata", "Aedes_aegypti.partial.gtf",
                    package="GenomicFeatures")
feature.type <- c("gene", "transcript", "exon", "CDS")
gr <- import(file, format="gtf", feature.type=feature.type)
txdb <- makeTxDbFromGRanges(gr)
txdb

}

