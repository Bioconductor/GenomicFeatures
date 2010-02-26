test_transcriptsBy <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                     package="GenomicFeatures"))

    checkException(transcriptsBy(data.frame()), silent = TRUE)
    checkException(transcriptsBy(txdb, "bad"), silent = TRUE)
    checkException(transcriptsBy(txdb, "tx"), silent = TRUE)

    dupCount <- function(x) {
        sum(sapply(x, function(y) anyDuplicated(values(y)[,"tx_id"])))
    }

    ## transcripts by gene
    txByGene <- transcriptsBy(txdb, "gene")
    checkTrue(validObject(txByGene))
    checkIdentical(dupCount(txByGene), 0L)
    checkIdentical(txByGene[[1]],
                   GRanges(seqnames = "chr21_random",
                           ranges   = IRanges(start=103280, end=164670),
                           strand   = strand("-"),
                           tx_name  = "uc002zka.1",
                           tx_id    = 120L,
                           gene_id  = "100132288"))

    ## transcripts by exon
    txByExon <- transcriptsBy(txdb, "exon")
    checkTrue(validObject(txByExon))
    checkIdentical(dupCount(txByExon), 0L)
    checkIdentical(txByExon[[1]],
                   GRanges(seqnames = "chr1",
                           ranges   = IRanges(start=4269, end=6628),
                           strand   = strand("-"),
                           tx_name  = "uc009vis.1",
                           tx_id    = 3L,
                           gene_id  = "653635"))

    ## transcripts by cds
    txByCds <- transcriptsBy(txdb, "cds")
    checkTrue(validObject(txByCds))
    checkIdentical(dupCount(txByCds), 0L)
    checkIdentical(txByCds[[1]],
                   GRanges(seqnames = "chr2",
                           ranges   = IRanges(start=31608, end=36385),
                           strand   = strand("-"),
                           tx_name  = "uc002qvt.1",
                           tx_id    = 4L,
                           gene_id  = "642273"))
}

test_exonsBy <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                     package="GenomicFeatures"))

    checkException(exonsBy(data.frame()), silent = TRUE)
    checkException(exonsBy(txdb, "bad"), silent = TRUE)
    checkException(exonsBy(txdb, "exon"), silent = TRUE)
    checkException(exonsBy(txdb, "cds"), silent = TRUE)

    dupCount <- function(x) {
        sum(sapply(x, function(y) anyDuplicated(values(y)[,"exon_id"])))
    }

    ## exons by transcript
    exonByTx <- exonsBy(txdb, "tx")
    checkTrue(validObject(exonByTx))
    checkIdentical(dupCount(exonByTx), 0L)
    checkIdentical(exonByTx[[2]],
                   GRanges(seqnames = c("chr1","chr1"),
                           ranges = IRanges(start = c(1116,2476),
                                              end = c(2090,4272)),
                           strand = strand(c("+","+")),
                           exon_name = as.character(c(NA,NA)),
                           exon_id = c(1L,4L) ))

    ## exons by gene
    exonByGene <- exonsBy(txdb, "gene")
    checkTrue(validObject(exonByGene))
    checkIdentical(dupCount(exonByGene), 0L)
}

test_cdsBy <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                     package="GenomicFeatures"))

    checkException(cdsBy(data.frame()), silent = TRUE)
    checkException(cdsBy(txdb, "bad"), silent = TRUE)
    checkException(cdsBy(txdb, "exon"), silent = TRUE)
    checkException(cdsBy(txdb, "cds"), silent = TRUE)

    dupCount <- function(x) {
        sum(sapply(x, function(y) anyDuplicated(values(y)[,"cds_id"])))
    }

    ## cds by transcript
    cdsByTx <- cdsBy(txdb, "tx")
    checkTrue(validObject(cdsByTx))
    checkIdentical(dupCount(cdsByTx), 0L)

    ## cds by gene
    cdsByGene <- cdsBy(txdb, "gene")
    checkTrue(validObject(cdsByGene))
    checkIdentical(dupCount(cdsByGene), 0L)
    checkIdentical(cdsByGene[[6]],
                   GRanges(seqnames = c("chr5","chr5"),
                           ranges = IRanges(start = c(258412,269844),
                                              end = c(259073,269964)),
                           strand = strand(c("-","-")),
                           cds_name = as.character(c(NA,NA)),
                           cds_id = c(53L,54L) ))
}
