test_transcriptsBy <- function()
{
    ## AN UNREALISTIC EDGE CASE
    ## ------------------------
    transcripts0 <- data.frame(
                        tx_id=c(26L, 5L, 11L),
                        tx_chrom=c("chr1", "chr2", "chr2"),
                        tx_strand=c("+", "-", "-"),
                        tx_start=c(1L, 16844685L, 16844685L),
                        tx_end=c(100L, 16844760L, 16844760L))
    splicings0 <-  data.frame(
                        tx_id=c(26L, 5L, 26L, 11L),
                        exon_rank=c(2L, 1L, 1L, 1L),
                        exon_start=c(1L, 16844685L, 1L, 16844685L),
                        exon_end=c(100L, 16844760L, 100L, 16844760L))
    txdb0 <- makeTranscriptDb(transcripts0, splicings0)

    ans <- transcriptsBy(txdb0, "exon")
    grg1 <- GRanges(seqnames="chr1",
                    ranges=IRanges(start=1, end=100),
                    strand=strand("+"),
                    tx_name=NA_character_,
                    tx_id=26L)
    grg2 <- GRanges(seqnames=c("chr2", "chr2"),
                    ranges=IRanges(start=c(16844685, 16844685),
                                   end=c(16844760,16844760)),
                    strand=strand(c("-", "-")),
                    tx_name=as.character(c(NA, NA)),
                    tx_id=c(5L, 11L))
    checkIdentical(ans, GRangesList(`1`=grg1, `2`=grg2))

    ans <- exonsBy(txdb0, "tx")
    grg5 <- GRanges(seqnames="chr2",
                    ranges=IRanges(start=16844685, end=16844760),
                    strand=strand("-"),
                    exon_name=NA_character_,
                    exon_id=2L)
    grg26 <- GRanges(seqnames="chr1",
                    ranges=IRanges(start=1, end=100),
                    strand=strand("+"),
                    exon_name=NA_character_,
                    exon_id=1L)[c(1L, 1L)]
    checkIdentical(ans, GRangesList(`5`=grg5, `11`=grg5, `26`=grg26))
                    
    ## WITH REAL DATA
    ## --------------
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
                           tx_id    = 120L))

    ## transcripts by exon
    txByExon <- transcriptsBy(txdb, "exon")
    checkTrue(validObject(txByExon))
    checkIdentical(dupCount(txByExon), 0L)
    checkIdentical(txByExon[[1]],
                   GRanges(seqnames = c("chr1", "chr1"),
                           ranges   = IRanges(start=1116, end=c(4121, 4272)),
                           strand   = strand(c("+", "+")),
                           tx_name  = c("uc001aaa.2", "uc009vip.1"),
                           tx_id    = c(1L, 2L)))

    ## transcripts by cds
    txByCds <- transcriptsBy(txdb, "cds")
    checkTrue(validObject(txByCds))
    checkIdentical(dupCount(txByCds), 0L)
    checkIdentical(txByCds[[1]],
                   GRanges(seqnames = "chr2",
                           ranges   = IRanges(start=31608, end=36385),
                           strand   = strand("-"),
                           tx_name  = "uc002qvt.1",
                           tx_id    = 4L))
}

test_exonsBy <- function()
{
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
