###

test_transcripts <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))

    checkException(transcripts(data.frame()), silent=TRUE)
    checkException(transcripts(txdb, filter=list(bad=1:10)), silent=TRUE)
    checkException(transcripts(txdb, columns="bad"), silent=TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    get_grg <- transcripts(txdb, filter=list(gene_id="139231"))
    want_grg <- GRanges(seqnames=factor("chrX", levels=seqlevels),
                        ranges=IRanges(start=c(103411156, 103430747),
                                       end  =c(103440582, 103440582)),
                        strand=strand("+"),
                        tx_id=142:143,
                        tx_name=c("uc004elw.3", "uc004elx.3"))
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

    filter <- list(tx_chrom=c("chr12", "chr14"), tx_strand="-")
    get_grg <- transcripts(txdb, columns=c("tx_id", "tx_name",
                                           "exon_id", "exon_rank"),
                                 filter=filter)
    want_grg <- GRanges(seqnames=factor("chr12", levels=seqlevels),
                        ranges=IRanges(start=52753790, end=52761309),
                        strand=strand("-"),
                        tx_id=87L,
                        tx_name="uc001sag.3",
                        exon_id=IntegerList(334:326),
                        exon_rank=IntegerList(1:9))
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

    filter <- list(gene_id=c("220004", "1183", "10186"))
    get_grg <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"),
                                 filter=filter)
    tx_id_col <- c(91L, 136:137)
    gene_id_col <- CharacterList("10186", "1183", "1183")
    want_grg <- GRanges(seqnames=factor(c("chr13", "chrX", "chrX"),
                                          levels=seqlevels),
                        ranges=IRanges(
                                   start=c(39917029, 10124985, 10124985),
                                   end  =c(40177356, 10205699, 10205699)),
                        strand=strand(c("-", "+", "+")),
                        tx_id=tx_id_col,
                        tx_name=c("uc001uxf.3", "uc004csy.4", "uc011mid.3"),
                        gene_id=gene_id_col)
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
    checkTrue(seqlevelsStyle(txdb)[1] == "NCBI")
    filter <- list(gene_id=c("220004", "1183", "10186"))
    get_grg <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"),
                                 filter=filter)
    checkTrue(all(seqnames(get_grg) == c("13", "X", "X")))

    tx_id_col <- c(91L, 136:137)
    gene_id_col <- CharacterList("10186", "1183", "1183")
    want_grg <- GRanges(seqnames=factor(c("13", "X", "X"),
                                          levels=seqlevels(txdb)),
                        ranges=IRanges(
                                   start=c(39917029, 10124985, 10124985),
                                   end  =c(40177356, 10205699, 10205699)),
                        strand=strand(c("-", "+", "+")),
                        tx_id=tx_id_col,
                        tx_name=c("uc001uxf.3", "uc004csy.4", "uc011mid.3"),
                        gene_id=gene_id_col)
    seqinfo(want_grg) <- seqinfo(txdb)
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)    
}

test_exons <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))

    checkException(exons(data.frame()), silent=TRUE)
    checkException(exons(txdb, filter=list(bad=1:10)), silent=TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    get_grg <- exons(txdb, filter=list(tx_name="uc001gde.2"))
    want_grg <- GRanges(seqnames=factor("chr1", levels=seqlevels),
                        ranges=IRanges(start=c(165513478, 165532742),
                                       end  =c(165514155, 165533185)),
                        strand=strand("+"),
                        exon_id=29:30)
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

    filter <- list(exon_chrom=c("chr5", "chr14"), exon_strand="-")
    get_grg <- exons(txdb, columns=c("exon_id", "tx_name", "gene_id"),
                           filter=filter)
    want_ranges <- IRanges(start=c(134363424, 134366966, 134369403, 170732985),
                           end  =c(134365011, 134367198, 134369964, 170735759))
    want_grg <- GRanges(seqnames=factor("chr5", levels=seqlevels),
                        ranges=want_ranges,
                        strand=strand("-"),
                        exon_id=182:185,
                        tx_name=CharacterList("uc010jea.3", "uc010jea.3",
                                              "uc010jea.3", "uc003mbe.2"),
                        gene_id=CharacterList("5307", "5307", "5307", NULL))
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)
}

test_cds <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                         package="GenomicFeatures"))

    checkException(cds(data.frame()), silent=TRUE)
    checkException(cds(txdb, filter=list(bad=1:10)), silent=TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    get_grg <- cds(txdb, filter=list(tx_name="uc001gde.2"))
    want_grg <- GRanges(seqnames=factor("chr1", levels=seqlevels),
                        ranges=IRanges(start=c(165513534, 165532742),
                                       end  =c(165514155, 165533061)),
                        strand=strand("+"),
                        cds_id=23:24)
    seqinfo(want_grg) <- seqinfo
    want_grg <- GenomicFeatures:::.assignMetadataList(want_grg, txdb)
    checkIdentical(get_grg, want_grg)

    filter <- list(cds_chrom=c("chr5", "chr14"), cds_strand="-")
    get_grg <- cds(txdb, columns=c("exon_id", "tx_name", "gene_id"),
                         filter=filter)
    want_ranges <- IRanges(start=c(134364469, 134366966, 134369403, 170735359),
                           end  =c(134365011, 134367198, 134369571, 170735634))
    want_grg <- GRanges(seqnames=factor("chr5", levels=seqlevels),
                        ranges=want_ranges,
                        strand=strand("-"),
                        exon_id=as(182:185, "IntegerList"),
                        tx_name=CharacterList("uc010jea.3", "uc010jea.3",
                                              "uc010jea.3", "uc003mbe.2"),
                        gene_id=CharacterList("5307", "5307", "5307", NULL))
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
    tx1 <- transcripts(txdb, columns=c("tx_id", "tx_name", "cds_id"))
    checkEquals(colnames(mcols(tx1)), c("tx_id", "tx_name", "cds_id"))
    tx2 <- transcripts(txdb, columns=c("TXID", "TXNAME", "CDSID"))
    checkEquals(colnames(mcols(tx2)), c("TXID", "TXNAME", "CDSID"))
    tx3 <- transcripts(txdb, columns=c(bob="CDSID"))
    checkEquals(colnames(mcols(tx3)), c("bob"))
    tx4 <- transcripts(txdb, columns=c(bob="cds_id"))
    checkEquals(colnames(mcols(tx4)), c("bob"))
    ## And these two cases should both explode. ;)
    checkException(transcripts(txdb, columns=c("")))
    checkException(transcripts(txdb, columns=c("bob")))
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
