###

test_transcripts <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))

    checkException(transcripts(data.frame()), silent = TRUE)
    checkException(transcripts(txdb, vals = list("bad" = 1:10)), silent = TRUE)
    checkException(transcripts(txdb, columns = "bad"), silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    get_grg <- transcripts(txdb, list("tx_id" = 3))
    want_grg <- GRanges(seqnames = factor("chr1", levels = seqlevels),
                        ranges = IRanges(start=145293233, end=145311154),
                        strand = strand("+"),
                        tx_id = 3L,
                        tx_name = "uc001emq.1")
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

    vals <- list(tx_chrom = c("chr11", "chr15"), tx_strand = "-")
    get_grg <- transcripts(txdb, vals, columns = c("tx_id","tx_name","exon_id"))
    want_ranges <- IRanges(start = c(53907623, 83509838),
                           end   = c(54051859, 83513971))
    want_grg <- GRanges(seqnames = factor("chr15", levels = seqlevels),
                        ranges = want_ranges,
                        strand = strand("-"),
                        tx_id = 71:72,
                        tx_name = c("uc010bfi.1", "uc002bjf.1"),
                        exon_id = IntegerList(617:633, 634))
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

    get_grg <- transcripts(txdb, vals=list(gene_id=c("220004", "1183", "10186")),
                                 columns=c("tx_id", "tx_name", "gene_id"))
    tx_id_col <- c(56L, 66L, 87L)
    gene_id_col <- CharacterList("220004", "10186", "1183")
    want_grg <- GRanges(seqnames = factor(c("chr11", "chr13", "chrX"),
                                          levels = seqlevels),
                        ranges = IRanges(
                                   start=c(61248585, 39917029, 10124985),
                                   end=c(61258400, 40177356, 10205699)),
                        strand = strand(c("+", "-", "+")),
                        tx_id = tx_id_col,
                        tx_name = c("uc001nru.2", "uc001uxf.3", "uc004csy.4"),
                        gene_id = gene_id_col)
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

}

## make sure that seqlevelsStyle behaves correctly
test_transcripts_seqlevelsStyleSwap <- function(){
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                               package="GenomicFeatures"))
    ## extra test to make sure that we can change names via seqlevelsStyle
    checkTrue(seqlevelsStyle(txdb) == "UCSC")
    seqlevelsStyle(txdb) <- "NCBI"
    checkTrue(seqlevelsStyle(txdb) == "NCBI")    
    get_grg <- transcripts(txdb, vals=list(gene_id=c("220004", "1183", "10186")),
                                 columns=c("tx_id", "tx_name", "gene_id"))
    checkIdentical(as.character(seqnames(get_grg)), c("11", "13", "X"))

    tx_id_col <- c(56L, 66L, 87L)
    gene_id_col <- CharacterList("220004", "10186", "1183")
    want_grg <- GRanges(seqnames = factor(c("11", "13", "X"),
                                          levels = seqlevels(txdb)),
                        ranges = IRanges(
                                   start=c(61248585, 39917029, 10124985),
                                   end=c(61258400, 40177356, 10205699)),
                        strand = strand(c("+", "-", "+")),
                        tx_id = tx_id_col,
                        tx_name = c("uc001nru.2", "uc001uxf.3", "uc004csy.4"),
                        gene_id = gene_id_col)
    seqinfo(want_grg) <- seqinfo(txdb)
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)    
}

test_exons <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))

    checkException(exons(data.frame()), silent = TRUE)
    checkException(exons(txdb, vals = list("bad" = 1:10)), silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    get_grg <- exons(txdb, list("exon_id" = 1))
    want_grg <- GRanges(seqnames = factor("chr1", levels = seqlevels),
                        ranges = IRanges(start=23853365, end=23855542),
                        strand = strand("+"),
                        exon_id = 1L)
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

    get_grg <- exons(txdb, vals = list(exon_chrom = c("chr22", "chr6_mcf_hap5"),
                                       exon_strand = "+"))
    want_ranges <- IRanges(start = c(3525991, 3527145, 3527255,
                                     3527481, 3527613, 3527833),
                           end   = c(3526257, 3527163, 3527367,
                                     3527535, 3527730, 3528397))
    want_grg <- GRanges(seqnames = factor("chr6_mcf_hap5", levels = seqlevels),
                        ranges = want_ranges,
                        strand = strand("+"),
                        exon_id = 864:869)
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)
}

test_cds <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                         package="GenomicFeatures"))

    checkException(cds(data.frame()), silent = TRUE)
    checkException(cds(txdb, vals = list("bad" = 1:10)), silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    get_grg <- cds(txdb, list(cds_id = 91))
    want_grg <- GRanges(seqnames = factor("chr2", levels = seqlevels),
                        ranges = IRanges(start=152285298, end=152285436),
                        strand = strand("+"),
                        cds_id = 91L)
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

    get_grg <- cds(txdb, vals = list(cds_chrom = c("chr22", "chr6_mcf_hap5"),
                                     cds_strand = "+"))
    want_grg <- GRanges(seqnames = factor("chr6_mcf_hap5", levels = seqlevels),
                        ranges = IRanges(start=c(3526118, 3527145, 3527255,
                                                 3527481, 3527613, 3527833),
                                         end=c(3526257, 3527163, 3527367,
                                               3527535, 3527730, 3527930)),
                        strand = strand("+"),
                        cds_id = 725:730)
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)
}

test_promoters <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                   package="GenomicFeatures"))
    tx <- transcripts(txdb)
    p <- suppressWarnings(promoters(txdb))
    checkEquals(c("tx_id", "tx_name"), colnames(mcols(p)))
    checkIdentical(promoters(tx[4:5]), p[4:5])
}


test_translateCols <- function(){
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
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
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                   package="GenomicFeatures"))
    de <- disjointExons(txdb, FALSE, TRUE)
    checkTrue(is(de$gene_id, "CharacterList"))
    checkTrue(is(de$tx_name, "CharacterList"))
    checkTrue(is(de$exonic_part, "integer"))
    checkIdentical(sum(elementLengths(de$gene_id) == 1), 864L)
    checkIdentical(sum(elementLengths(de$gene_id) == 2), 0L)

    de <- disjointExons(txdb, TRUE, FALSE)
    checkTrue(all(names(mcols(de)) %in% c("gene_id", "exonic_part")))
    checkIdentical(sum(elementLengths(de$gene_id) == 1), 864L)
    checkIdentical(sum(elementLengths(de$gene_id) == 2), 0L)
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
