### =========================================================================
### makeTxDbFromGRanges()
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the 'transcripts' data frame
###

.extract_transcripts_from_GRanges <- function(gr, type, ID, Name)
{
    idx <- which(type == "mRNA")
    tx_id <- ID[idx]
    if (anyDuplicated(tx_id))
        stop("transcript IDs are not unique")
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
    tmp <- IntegerList(unname(split(rank, parent_id)))
    if (any(any(duplicated(tmp))))
        return(NULL)
    if (for.cds)
        return(rank)
    if (!(all(min(tmp) == 1L) && all(max(tmp) == elementLengths(tmp))))
        return(NULL)
    rank
}

### Can be used to extract exons or cds.
.extract_exons <- function(gr, type, ID, Name, Parent, for.cds=FALSE)
{
    feature.type <- if (for.cds) "CDS" else "exon"
    exon_idx <- which(type == feature.type)
    exon_Parent <- Parent[exon_idx]
    tx_id <- factor(unlist(exon_Parent, use.names=FALSE))
    nparent_per_ex <- elementLengths(exon_Parent)
    exon_id <- rep.int(ID[exon_idx], nparent_per_ex)

    exon_rank <- .extract_rank_from_id(exon_id, tx_id, for.cds=for.cds)
    if (is.null(exon_rank))
        stop("failed to extract the ", feature.type, " ranks ",
             "from the ", feature.type, " IDs")

    exon_name <- rep.int(Name[exon_idx], nparent_per_ex)
    exon_chrom <- rep.int(seqnames(gr)[exon_idx], nparent_per_ex)
    exon_strand <- rep.int(strand(gr)[exon_idx], nparent_per_ex)
    exon_start <- rep.int(start(gr)[exon_idx], nparent_per_ex)
    exon_end <- rep.int(end(gr)[exon_idx], nparent_per_ex)
    ans <- data.frame(
        tx_id=tx_id,
        exon_rank=exon_rank,
        exon_name=exon_name,
        exon_chrom=exon_chrom,
        exon_strand=exon_strand,
        exon_start=exon_start,
        exon_end=exon_end,
        stringsAsFactors=FALSE
    )
    if (for.cds)
        colnames(ans) <- sub("^exon", "cds", colnames(ans))
    ans
}

.map_cds_to_exon <- function(cds, exons)
{
    query <- GRanges(cds$tx_id, IRanges(cds$cds_start, cds$cds_end))
    subject <- GRanges(exons$tx_id, IRanges(exons$exon_start, exons$exon_end))
    cds2exon <- findOverlaps(query, subject, type="within")
    if (anyDuplicated(queryHits(cds2exon)))
        stop("some exons contain more than one CDS")
    selectHits(cds2exon, select="arbitrary")
}

.extract_splicings_from_GRanges <- function(gr, type, ID, Name, Parent)
{
    exons <- .extract_exons(gr, type, ID, Name, Parent)
    cds <- .extract_exons(gr, type, ID, Name, Parent, for.cds=TRUE)

    cds_tx_id <- factor(cds$tx_id, levels=levels(exons$tx_id))
    if (any(is.na(cds_tx_id)))
        stop("some CDS cannot be mapped to an exon")
    cds$tx_id <- cds_tx_id

    cds2exon <- .map_cds_to_exon(cds, exons)
    if (any(is.na(cds2exon)))
        stop("some CDS cannot be mapped to an exon")
    if (any(cds$cds_chrom != exons$exon_chrom[cds2exon])
     || any(cds$cds_strand != exons$exon_strand[cds2exon]))
        stop("some CDS are on a chrom/strand that differs from ",
             "the chrom/strand of the exon they're mapped to")
    cds_name <- rep.int(NA_character_, nrow(exons))
    cds_start <- cds_end <- rep.int(NA_integer_, nrow(exons))
    cds_name[cds2exon] <- cds$cds_name
    cds_start[cds2exon] <- cds$cds_start
    cds_end[cds2exon] <- cds$cds_end
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
        stop("some exons are mapped to transcripts that cannot be found")
    splicings$tx_id <- splicings_tx_id
    genes_tx_id <- match(genes$tx_id, transcripts$tx_id)
    if (any(is.na(genes_tx_id)))
        stop("some genes are mapped to transcripts that cannot be found")
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

