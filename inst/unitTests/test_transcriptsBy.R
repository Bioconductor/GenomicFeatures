test_transcriptsBy <- function()
{
    ## A TOY CASE
    ## ----------
    transcripts0 <- data.frame(
                        tx_id=c(26L, 5L, 11L),
                        tx_name=c("A", "B", "C"),
                        tx_chrom=c("chr1", "chr2", "chr2"),
                        tx_strand=c("+", "-", "-"),
                        tx_start=c(1L, 16844685L, 16844685L),
                        tx_end=c(100L, 16844760L, 16844760L))
    splicings0 <-  data.frame(
                        tx_id=c(26L, 5L, 26L, 11L),
                        exon_rank=c(2L, 1L, 1L, 1L),
                        exon_start=c(1L, 16844685L, 1L, 16844685L),
                        exon_end=c(100L, 16844760L, 100L, 16844760L))

    suppressWarnings(txdb0 <- makeTxDb(transcripts0, splicings0))

    seqinfo <- seqinfo(txdb0)
    seqlevels <- seqlevels(seqinfo)

    ans <- transcriptsBy(txdb0, by="exon")
    grg1 <- GRanges(seqnames=factor("chr1", levels = seqlevels),
                    ranges=IRanges(start=c(1, 1), end=100),
                    strand=strand("+"),
                    tx_id=c(26L, 26L),
                    tx_name="A",
                    exon_rank=1:2)
    grg2 <- GRanges(seqnames=factor(c("chr2", "chr2"), levels = seqlevels),
                    ranges=IRanges(start=c(16844685, 16844685),
                                   end=c(16844760,16844760)),
                    strand=strand(c("-", "-")),
                    tx_id=c(5L, 11L),
                    tx_name=c("B", "C"),
                    exon_rank=1L)
    want <- GRangesList(`1`=grg1, `2`=grg2)
    seqinfo(want) <- seqinfo
    want <- GenomicFeatures:::.assignMetadataList(want, txdb0)
    checkIdentical(ans, want)

    ans <- exonsBy(txdb0, "tx")
    grg5 <- GRanges(seqnames=factor("chr2", levels=seqlevels),
                    ranges=IRanges(start=16844685, end=16844760),
                    strand=strand("-"),
                    exon_id=2L,
                    exon_name=NA_character_,
                    exon_rank=1L)
    grg26 <- GRanges(seqnames=factor("chr1", levels=seqlevels)[c(1L, 1L)],
                    ranges=IRanges(start=1, end=100),
                    strand=strand("+"),
                    exon_id=1L,
                    exon_name=NA_character_,
                    exon_rank=1:2)
    want <- GRangesList(`5`=grg5, `11`=grg5, `26`=grg26)
    seqinfo(want) <- seqinfo
    want <- GenomicFeatures:::.assignMetadataList(want, txdb0)
    checkIdentical(ans, want)
                    
    ## WITH REAL DATA
    ## --------------
    txdb1 <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite",
                                     package="GenomicFeatures"))

    checkException(transcriptsBy(data.frame()), silent = TRUE)
    checkException(transcriptsBy(txdb1, by="bad"), silent = TRUE)
    checkException(transcriptsBy(txdb1, by="tx"), silent = TRUE)

    seqinfo <- seqinfo(txdb1)
    want <- GenomicFeatures:::.assignMetadataList(want, txdb1)
    seqlevels <- seqlevels(seqinfo)

    dupCount <- function(x) {
        sum(sapply(x, function(y) anyDuplicated(mcols(y)[,"tx_id"])))
    }

    ## transcripts by gene
    txByGene <- transcriptsBy(txdb1, by="gene")
    checkTrue(validObject(txByGene))
    checkIdentical(dupCount(txByGene), 0L)
    want <- GRanges(seqnames = factor("chr6", levels=seqlevels),
                    ranges   = IRanges(start=c(10412551, 10414300, 10414300),
                                       end  =c(10416402, 10416190, 10416402)),
                    strand   = strand("+"),
                    tx_id    = 46:48,
                    tx_name  = c("uc003myw.3", "uc003myy.1", "uc003myx.3"))
    seqinfo(want) <- seqinfo
#    metadata(want)[[1]] <- DataFrame(metadata(txdb1)) ##WTH?
    checkIdentical(txByGene[["100130275"]], want)

    ## transcripts by exon
    txByExon <- transcriptsBy(txdb1, by="exon")
    checkTrue(validObject(txByExon))
    checkIdentical(dupCount(txByExon), 0L)
    want <- GRanges(seqnames  = factor("chr1", levels=seqlevels),
                    ranges    = IRanges(start=c(32671236, 32671236),
                                        end  =c(32674288, 32674288)),
                    strand    = strand("+"),
                    tx_id     = c(1L, 3L),
                    tx_name   = c("uc001bum.2", "uc010ogz.1"),
                    exon_rank = 1L)
    seqinfo(want) <- seqinfo
#    metadata(want)[[1]] <- DataFrame(metadata(txdb1)) ##WTH?
    checkIdentical(txByExon[[1L]], want)

    ## transcripts by cds
    txByCds <- transcriptsBy(txdb1, by="cds")
    checkTrue(validObject(txByCds))
    checkIdentical(dupCount(txByCds), 0L)
    want <- GRanges(seqnames  = factor("chr1", levels=seqlevels),
                    ranges    = IRanges(start=32671236, end=32674288),
                    strand    = strand("+"),
                    tx_id     = 1L,
                    tx_name   = "uc001bum.2",
                    exon_rank = 1L)
    seqinfo(want) <- seqinfo
#    metadata(want)[[1]] <- DataFrame(metadata(txdb1)) ##WTH?
    checkIdentical(txByCds[[1L]], want)

    ## threeUTRsByTranscript, fiveUTRsByTranscript
    threeUTRs <- threeUTRsByTranscript(txdb1)
    checkTrue(validObject(threeUTRs))
    fiveUTRs <- fiveUTRsByTranscript(txdb1)
    checkTrue(validObject(fiveUTRs))
}

## make sure that seqlevelsStyle behaves correctly
test_transcriptsBy_seqlevelsStyleSwap <- function(){
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                               package="GenomicFeatures"))    
    get_grl <- transcriptsBy(txdb, by="gene")
    checkTrue(seqlevelsStyle(txdb) == "UCSC")
    checkTrue(all(seqnames(get_grl[["100130275"]]) =="chr6"))
    
    seqlevelsStyle(txdb) <- "NCBI"
    checkTrue(seqlevelsStyle(txdb)[1] == "NCBI")    
    get_grlN <- transcriptsBy(txdb, by="gene")
    checkTrue(all(seqnames(get_grlN[["100130275"]]) == "6"))
}

test_exonsBy <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite",
                                     package="GenomicFeatures"))

    checkException(exonsBy(data.frame()), silent = TRUE)
    checkException(exonsBy(txdb, "bad"), silent = TRUE)
    checkException(exonsBy(txdb, "exon"), silent = TRUE)
    checkException(exonsBy(txdb, "cds"), silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    dupCount <- function(x) {
        sum(sapply(x, function(y) anyDuplicated(mcols(y)[,"exon_id"])))
    }

    ## exons by transcript
    exonByTx <- exonsBy(txdb, "tx")
    checkTrue(validObject(exonByTx))
    checkIdentical(dupCount(exonByTx), 0L)
    want <- GRanges(seqnames = factor("chr1", levels=seqlevels),
                    ranges = IRanges(start = c(153330330, 153330745, 153333120),
                                     end = c(153330357, 153330909, 153333503)),
                    strand = strand(c("+", "+", "+")),
                    exon_id = 8:10,
                    exon_name = as.character(c(NA, NA, NA)),
                    exon_rank = 1:3)
    seqinfo(want) <- seqinfo
#    metadata(want)[[1]] <- DataFrame(metadata(txdb))  ##WTH?
    checkIdentical(exonByTx[[4L]], want)

    ## exons by gene
    exonByGene <- exonsBy(txdb, "gene")
    checkTrue(validObject(exonByGene))
    checkIdentical(dupCount(exonByGene), 0L)
}

test_cdsBy <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite",
                                     package="GenomicFeatures"))

    checkException(cdsBy(data.frame()), silent = TRUE)
    checkException(cdsBy(txdb, "bad"), silent = TRUE)
    checkException(cdsBy(txdb, "exon"), silent = TRUE)
    checkException(cdsBy(txdb, "cds"), silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    dupCount <- function(x) {
        sum(sapply(x, function(y) anyDuplicated(mcols(y)[,"cds_id"])))
    }

    ## cds by transcript
    cdsByTx <- cdsBy(txdb, "tx")
    checkTrue(validObject(cdsByTx))
    checkIdentical(dupCount(cdsByTx), 0L)

    ## cds by gene
    cdsByGene <- cdsBy(txdb, "gene")
    checkTrue(validObject(cdsByGene))
    checkIdentical(dupCount(cdsByGene), 0L)
    want <- GRanges(seqnames = factor("chr13", levels = seqlevels),
                    ranges = IRanges(start = c(39918073, 39952565, 40174969),
                                     end = c(39918191, 39952663, 40175353)),
                    strand = strand("-"),
                    cds_id = 279:281,
                    cds_name = as.character(c(NA, NA, NA)))
    seqinfo(want) <- seqinfo
#    metadata(want)[[1]] <- DataFrame(metadata(txdb))  ##WTH?
    checkIdentical(cdsByGene[["10186"]], want)
}

test_intronsByTranscript <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite",
                                     package="GenomicFeatures"))

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    intronByTx <- intronsByTranscript(txdb)
    checkTrue(validObject(intronByTx))
    want <- GRanges(seqnames = factor("chr1", levels = seqlevels),
                    ranges = IRanges(start=c(153330358, 153330910),
                                     end=c(153330744, 153333119)),
                    strand = strand("+"))
    seqinfo(want) <- seqinfo
    checkIdentical(intronByTx[[4L]], want)
}
