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

## make sure that seqnameStyle behaves correctly
test_transcripts_seqnameStyleSwap <- function(){
    txdb <- loadDb(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))
    ## extra test to make sure that we can change names via seqnameStyle
    checkTrue(seqnameStyle(txdb) == "UCSC")
    seqnameStyle(txdb) <- "NCBI"
    checkTrue(seqnameStyle(txdb) == "NCBI")    
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

