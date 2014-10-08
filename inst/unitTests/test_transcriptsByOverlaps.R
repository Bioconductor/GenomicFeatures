###

test_transcriptsByOverlaps <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))

    checkException(transcriptsByOverlaps(txdb), silent = TRUE)
    checkException(transcriptsByOverlaps(txdb, IRanges()), silent = TRUE)
    checkException(transcriptsByOverlaps(txdb, GRanges(), columns = "bad"),
                   silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    gr <- GRanges(seqnames = "chr8",
                  ranges = IRanges(start=1, end=29100000),
                  strand = strand("-"))
    want <-
      GRanges(seqnames = factor("chr8", levels = seqlevels),
              ranges = IRanges(start=29003375, end=29120610),
              strand = strand("-"),
              tx_id = 43L,
              tx_name = "uc003xhj.2")
    seqinfo(want) <- seqinfo
    want <- GenomicFeatures:::.assignMetadataList(want, txdb)
    checkIdentical(transcriptsByOverlaps(txdb, gr), want)

    ranges <- IRanges(start = c(1, 150000000, 1, 150000000),
                      end   = c(35000000, 250000000, 35000000, 250000000))
    chrom <- c("chr1", "chr1", "chr2", "chr2")
    strand <- strand(c("+", "-", "+", "-"))
    gr <- GRanges(seqnames = chrom, ranges = ranges, strand = strand)
    want <- GRanges(seqnames =
                    factor(c("chr1", "chr1", "chr1", "chr2", "chr2"),
                           levels = seqlevels),
                    ranges = IRanges(start = c(23853365, 32671236, 150547027,
                                               10281981, 166626170),
                                     end   = c(23855542, 32674288, 150552214,
                                               10350942, 166650803)),
                    strand = strand(c("+", "+", "-", "+", "-")),
                    tx_id = c(1L, 2L, 11L, 12L, 17L))
    seqinfo(want) <- seqinfo
    want <- GenomicFeatures:::.assignMetadataList(want, txdb)
    checkIdentical(transcriptsByOverlaps(txdb, gr, columns="tx_id"), want)
}

test_exonsByOverlaps <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))

    checkException(exonsByOverlaps(txdb), silent = TRUE)
    checkException(exonsByOverlaps(txdb, IRanges()), silent = TRUE)
    checkException(exonsByOverlaps(txdb, GRanges(), columns = "bad"),
                   silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    gr <- GRanges(seqnames = "chr8",
                  ranges = IRanges(start=129000000, end=141000000),
                  strand = strand("-"))
    want <-
      GRanges(seqnames = factor(rep("chr8", 3), levels = seqlevels),
              ranges = IRanges(start=c(129022628, 140922391, 140998934),
                               end  =c(129022857, 140922544, 140999044)),
              strand = strand(rep("-", 3)),
              exon_id = 359:361)
    seqinfo(want) <- seqinfo
    want <- GenomicFeatures:::.assignMetadataList(want, txdb)
    checkIdentical(exonsByOverlaps(txdb, gr), want)
}

test_cdsByOverlaps <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))

    checkException(cdsByOverlaps(txdb), silent = TRUE)
    checkException(cdsByOverlaps(txdb, IRanges()), silent = TRUE)
    checkException(cdsByOverlaps(txdb, GRanges(), columns = "bad"),
                   silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    gr <- GRanges(seqnames = "chr8",
                  ranges = IRanges(start=129000000, end=141000000),
                  strand = strand("-"))
    want <-
      GRanges(seqnames = factor(rep("chr8", 2), levels = seqlevels),
              ranges = IRanges(start=c(140922391, 140998934),
                               end  =c(140922544, 140999044)),
              strand = strand(rep("-", 2)),
              cds_id = 304:305)
    seqinfo(want) <- seqinfo
    want <- GenomicFeatures:::.assignMetadataList(want, txdb)
    checkIdentical(cdsByOverlaps(txdb, gr), want)
}
