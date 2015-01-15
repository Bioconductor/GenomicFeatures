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
### Extract the 'transcripts' data frame
###

.extract_transcripts_from_GRanges <- function(gr, type, ID, Name)
{
    idx <- which(type == "mRNA")
    tx_id <- ID[idx]
    if (anyDuplicated(tx_id))
        stop(wmsg("transcript IDs are not unique"))
    data.frame(
        tx_id=tx_id,
        tx_name=Name[idx],
        tx_chrom=seqnames(gr)[idx],
        tx_strand=strand(gr)[idx],
        tx_start=start(gr)[idx],
        tx_end=end(gr)[idx],
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
    if (any(elementLengths(tx_chrom) > 1L))
        stop(wmsg("Some transcripts have exons on more than one chromosome. ",
                  "Cannot infer the exon ranks."))

    strand_per_tx <- split(Rle(exon_strand), tx_id)
    tx_strand <- runValue(strand_per_tx)
    if (any(elementLengths(tx_strand) > 1L))
        stop(wmsg("Some transcripts have exons on both strands. ",
                  "Cannot infer the exon ranks."))
    minus_idx <- which(as.character(tx_strand) == "-")

    ex_per_tx <- split(IRanges(exon_start, exon_end), tx_id)
    if (!identical(elementLengths(reduce(ex_per_tx)),
                   elementLengths(ex_per_tx)))
        stop(wmsg("Some transcripts have exons not separated by introns. ",
                  "Cannot infer the exon ranks."))

    start_per_tx <- start(ex_per_tx)
    start_per_tx[minus_idx] <- start_per_tx[minus_idx] * (-1L)
    rank <- .rank.CompressedList(start_per_tx, ties.method="first")
    unsplit(rank, tx_id)
}

### Can be used to extract exons or cds.
.extract_exons <- function(gr, type, ID, Name, Parent, for.cds=FALSE)
{
    feature.type <- if (for.cds) "CDS" else "exon"
    exon_idx <- which(type == feature.type)
    exon_Parent <- Parent[exon_idx]
    if (any(any(duplicated(exon_Parent))))
        stop(wmsg("some ", feature.type, "s are mapped twice ",
                  "to the same transcript"))

    tx_id <- factor(unlist(exon_Parent, use.names=FALSE))
    nparent_per_ex <- elementLengths(exon_Parent)
    exon_id <- rep.int(ID[exon_idx], nparent_per_ex)
    exon_name <- rep.int(Name[exon_idx], nparent_per_ex)
    exon_chrom <- rep.int(seqnames(gr)[exon_idx], nparent_per_ex)
    exon_strand <- rep.int(strand(gr)[exon_idx], nparent_per_ex)
    exon_start <- rep.int(start(gr)[exon_idx], nparent_per_ex)
    exon_end <- rep.int(end(gr)[exon_idx], nparent_per_ex)

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
        stop(wmsg("the following transcripts have CDS that are mapped ",
                  "to more than one exon: ", bad_tx))
    }

    cds2exon <- selectHits(hits, select="arbitrary")
    bad_tx <- unique(cds$tx_id[is.na(cds2exon)])
    if (length(bad_tx)) {
        bad_tx <- paste0(bad_tx, collapse=", ")
        stop(wmsg("the following transcripts have CDS that cannot ",
                  "be mapped to an exon: ", bad_tx))
    }

    bad_tx <- unique(exons$tx_id[s_hits[duplicated(s_hits)]])
    if (length(bad_tx)) {
        bad_tx <- paste0(bad_tx, collapse=", ")
        warning(wmsg("the following transcripts have exons containing ",
                     "more than one CDS (only the first CDS was kept ",
                     "for each exon): ", bad_tx))
    }
    selectHits(t(hits), select="first")
}

.extract_splicings_from_GRanges <- function(gr, type, ID, Name, Parent)
{
    exons <- .extract_exons(gr, type, ID, Name, Parent)
    cds <- .extract_exons(gr, type, ID, Name, Parent, for.cds=TRUE)

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

.extract_genes_from_GRanges <- function(type, ID, Name, Parent)
{
    idx <- which(type == "gene")
    geneID2Name <- Name[idx]
    names(geneID2Name) <- ID[idx]

    idx <- which(type == "mRNA")
    tx2genes <- Parent[idx]
    tx_id <- rep.int(factor(ID[idx]), elementLengths(tx2genes))
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

### Metadata columns we care about:
### o If GRanges obtained from GFF3: type, ID, Name, Parent,
### o If GRanges obtained from GTF:
###     type, gene_id, transcript_id, exon_number, gene_name, transcript_name
makeTxDbFromGRanges <- function(gr, metadata=NULL)
{
    gr_mcols <- mcols(gr)
    type <- gr_mcols[ , "type"]
    ID <- gr_mcols[ , "ID"]
    Name <- gr_mcols[ , "Name"]
    Parent <- gr_mcols[ , "Parent"]

    transcripts <- .extract_transcripts_from_GRanges(gr, type, ID, Name)
    splicings <- .extract_splicings_from_GRanges(gr, type, ID, Name, Parent)
    genes <- .extract_genes_from_GRanges(type, ID, Name, Parent)
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
library(rtracklayer)

## Test with GRanges obtained from GFF3 files
## ==========================================

feature.type <- c("gene", "mRNA", "exon", "CDS")

file <- system.file("extdata", "GFF3_files", "TheCanonicalGene_v1.gff3",
                    package="GenomicFeatures")
gr <- import(file, format="gff3", feature.type=feature.type)
txdb <- makeTxDbFromGRanges(gr)
txdb

file <- system.file("extdata", "a.gff3", package="GenomicFeatures")
gr <- import(file, format="gff3", feature.type=feature.type)
txdb <- makeTxDbFromGRanges(gr)
txdb

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

