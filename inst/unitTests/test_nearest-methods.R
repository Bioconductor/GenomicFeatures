quiet <- suppressWarnings

test_GenomicRanges_distance <- function()
{
    genes <- data.frame(
        tx_id=1:3,
        gene_id=c("gene1", "gene1", "gene2"))
    transcripts <- data.frame(
        tx_id=1:3,
        tx_chrom="chr1",
        tx_strand=c("+", "+", "-"),
        tx_start=c(1, 2001, 3001),
        tx_end=c(999, 2199, 3199))
    splicings <-  data.frame(
        tx_id=c(1L, 2L, 2L, 2L, 3L, 3L),
        cds_id=c(10L, 11L, 12L, 13L, 14L, 15L),
        exon_rank=c(1, 1, 2, 3, 1, 2),
        exon_start=c(1, 2001, 2101, 2131, 3001, 3131),
        exon_end=c(999, 2085, 2144, 2199, 3085, 3199),
        cds_start=c(1, 2022, 2101, 2131, 3001, 3131),
        cds_end=c(999, 2085, 2144, 2193, 3085, 3199))
    txdb <- quiet(makeTxDb(transcripts, splicings, genes))

    gr <- GRanges("chr1", IRanges(1050, width=1))
    strand(gr) <- "-"
    d <- quiet(distance(gr, txdb, id="gene1", type="gene"))
    checkTrue(is.na(d))
    strand(gr) <- "+"
    d_pos <- quiet(distance(gr, txdb, id="gene1", type="gene"))
    strand(gr) <- "*"
    d_star <- quiet(distance(gr, txdb, id="gene1", type="gene"))
    checkIdentical(d_pos, d_star)

    d_tx <- quiet(distance(gr, txdb, id="3", type="tx"))
    d_cds <- quiet(distance(gr, txdb, id="14", type="cds"))
    checkIdentical(d_tx, d_cds)
}
