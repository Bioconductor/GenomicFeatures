###

test_transcripts <- function()
{
    txdb <- loadDb(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))

    checkException(transcripts(data.frame()), silent = TRUE)
    checkException(transcripts(txdb, vals = list("bad" = 1:10)), silent = TRUE)
    checkException(transcripts(txdb, columns = "bad"), silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    get_grg <- transcripts(txdb, list("tx_id" = 3))
    want_grg <- GRanges(seqnames = factor("chr1", levels = seqlevels),
                        ranges = IRanges(start=4269, end=6628),
                        strand = strand("-"),
                        tx_id = 3L,
                        tx_name = "uc009vis.1")
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

    vals <- list(tx_chrom = c("chr1", "chr5"), tx_strand = "-")
    get_grg <- transcripts(txdb, vals, columns = c("tx_id","tx_name","exon_id"))
    want_ranges <- IRanges(start = c(4269, 257875),
                           end   = c(6628, 271297))
    want_grg <- GRanges(seqnames = factor(c("chr1", "chr5"), levels = seqlevels),
                        ranges = want_ranges,
                        strand = strand(rep("-", 2)),
                        tx_id = c(3L, 15L),
                        tx_name = c("uc009vis.1", "uc003jam.1"),
                        exon_id = IntegerList(c(5,6,7,8), c(77,78,79)))
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

    get_grg <- transcripts(txdb, vals=list(gene_id=c("3081", "9501")),
                                 columns=c("tx_id", "tx_name", "gene_id"))
    tx_id_col <- c(49:51, 83:84)
    gene_id_col <- CharacterList(as.list(c(rep.int("9501", 3),
                                           rep.int("3081", 2))))
    want_grg <- GRanges(seqnames = factor(c(rep.int("chr17", 3),
                                            rep.int("chr3_random", 2)),
                                          levels = seqlevels),
                        ranges = IRanges(
                                   start=c(rep.int(62294, 3), 18988, 49418),
                                   end=c(177378, 202576, 202576, 73308, 73308)),
                        strand = strand(c("-", "-", "-", "+", "+")),
                        tx_id = tx_id_col,
                        tx_name = c("uc002frd.1", "uc002fre.1", "uc002frf.1",
                                    "uc003fzi.1", "uc003fzj.1"),
                        gene_id = gene_id_col)
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

}

## make sure that seqlevelStyle behaves correctly
test_transcripts_seqlevelStyleSwap <- function(){
    txdb <- loadDb(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))
    ## extra test to make sure that we can change names via seqlevelStyle
    checkTrue(seqlevelStyle(txdb) == "UCSC")
    seqlevelStyle(txdb) <- "NCBI"
    checkTrue(seqlevelStyle(txdb) == "NCBI")    
    get_grg <- transcripts(txdb, vals=list(gene_id=c("3081", "9501")),
                                 columns=c("tx_id", "tx_name", "gene_id"))
    checkTrue(as.character(seqnames(get_grg))[1] =="17")
    tx_id_col <- c(49:51, 83:84)
    gene_id_col <- CharacterList(as.list(c(rep.int("9501", 3),
                                           rep.int("3081", 2))))
    want_grg <- GRanges(seqnames = factor(c(rep.int("17", 3),
                                            rep.int("chr3_random", 2)),
                                          levels = seqlevels(seqinfo(txdb))),
                        ranges = IRanges(
                                   start=c(rep.int(62294, 3), 18988, 49418),
                                   end=c(177378, 202576, 202576, 73308, 73308)),
                        strand = strand(c("-", "-", "-", "+", "+")),
                        tx_id = tx_id_col,
                        tx_name = c("uc002frd.1", "uc002fre.1", "uc002frf.1",
                                    "uc003fzi.1", "uc003fzj.1"),
                        gene_id = gene_id_col)
    seqinfo(want_grg) <- seqinfo(txdb)
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)    
}

test_exons <- function()
{
    txdb <- loadDb(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))

    checkException(exons(data.frame()), silent = TRUE)
    checkException(exons(txdb, vals = list("bad" = 1:10)), silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    get_grg <- exons(txdb, list("exon_id" = 1))
    want_grg <- GRanges(seqnames = factor("chr1", levels = seqlevels),
                        ranges = IRanges(start=1116, end=2090),
                        strand = strand("+"),
                        exon_id = 1L)
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

    get_grg <- exons(txdb, vals = list(exon_chrom = c("chr1", "chr5"),
                                       exon_strand = "-"))
    want_ranges <- IRanges(start = c(4269,4833,5659,6470,257875,269844,271208),
                           end   = c(4692,4901,5805,6628,259073,269974,271297))
    want_grg <- GRanges(seqnames =
                          factor(c(rep("chr1", 4), rep("chr5", 3)), levels = seqlevels),
                        ranges = want_ranges,
                        strand = strand(rep("-", 7)),
                        exon_id = c(5L, 6L, 7L, 8L, 77L, 78L, 79L))
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)
}

test_cds <- function()
{
    txdb <- loadDb(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                         package="GenomicFeatures"))

    checkException(cds(data.frame()), silent = TRUE)
    checkException(cds(txdb, vals = list("bad" = 1:10)), silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    get_grg <- cds(txdb, list("cds_id" = 91))
    want_grg <- GRanges(seqnames = factor("chr10", levels = seqlevels),
                        ranges = IRanges(start=82997, end=84054),
                        strand = strand("-"),
                        cds_id = 91L)
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

    get_grg <- cds(txdb, vals = list(cds_chrom = c("chr1", "chr5"),
                                     cds_strand = "-"))
    want_grg <- GRanges(seqnames = factor(c("chr5", "chr5"), levels = seqlevels),
                        ranges = IRanges(start=c(258412,269844), end=c(259073,269964)),
                        strand = strand(c("-", "-")),
                        cds_id = 53:54)
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)
}

test_promoters <- function()
{
    txdb <- loadDb(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                   package="GenomicFeatures"))
    tx <- transcripts(txdb)
    p <- suppressWarnings(promoters(txdb))
    checkEquals(c("tx_id", "tx_name"), colnames(mcols(p)))
    checkIdentical(promoters(tx[4:5]), p[4:5])
}


test_translateCols <- function(){
    txdb <- loadDb(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                   package="GenomicFeatures"))
    tx1 <- transcripts(txdb, columns = c("tx_id", "tx_name", "cds_id"))
    checkEquals(colnames(mcols(tx1)), c("tx_id", "tx_name", "cds_id"))
    tx2 <- transcripts(txdb, columns = c("TXID", "TXNAME", "CDSID"))
    checkEquals(colnames(mcols(tx2)), c("TXID", "TXNAME", "CDSID"))
    tx3 <- transcripts(txdb, columns = c(bob="CDSID"))
    checkEquals(colnames(mcols(tx3)), c("bob"))
    tx4 <- transcripts(txdb, columns = c(bob="cds_id"))
    checkEquals(colnames(mcols(tx4)), c("bob"))
    ## And these two cases should both explode. ;)
    checkException(transcripts(txdb, columns = c("")))
    checkException(transcripts(txdb, columns = c("bob")))
}

test_disjointExons <- function()
{
    txdb <- loadDb(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                   package="GenomicFeatures"))
    de <- disjointExons(txdb, FALSE, TRUE)
    checkTrue(is(de$gene_id, "CharacterList"))
    checkTrue(is(de$tx_name, "CharacterList"))
    checkTrue(is(de$exonic_part, "integer"))
    checkIdentical(sum(elementLengths(de$gene_id) == 1), 371L)
    checkIdentical(sum(elementLengths(de$gene_id) == 2), 0L)

    de <- disjointExons(txdb, TRUE, FALSE)
    checkTrue(all(names(mcols(de)) %in% c("gene_id", "exonic_part")))
    checkIdentical(sum(elementLengths(de$gene_id) == 1), 361L)
    checkIdentical(sum(elementLengths(de$gene_id) == 2), 18L)
}

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
    txdb <- quiet(makeTranscriptDb(transcripts, splicings, genes))

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
