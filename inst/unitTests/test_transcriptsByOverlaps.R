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

    gr <- GRanges(seqnames = "chrX",
                  ranges = IRanges(start=54071000, width=1),
                  strand = strand("-"))
    want <-
      GRanges(seqnames = factor("chrX", levels = seqlevels),
              ranges = IRanges(start=53963113, end=54071569),
              strand = strand("-"),
              tx_id = 147L,
              tx_name = "uc004dsu.3")
    seqinfo(want) <- seqinfo
    want <- GenomicFeatures:::.assignMetadataList(want, txdb)
    checkIdentical(transcriptsByOverlaps(txdb, gr), want)

    ranges <- IRanges(start = c(113000000, 54071000, 54071000),
                      width = c(  5000000,        1,        1))
    chrom <- c("chr3", "chrX", "chrX")
    strand <- strand(c("+", "+", "-"))
    gr <- GRanges(seqnames = chrom, ranges = ranges, strand = strand)
    want <- GRanges(seqnames =
                    factor(c("chr3", "chr3", "chr3", "chrX"),
                           levels = seqlevels),
                    ranges = IRanges(start = c(113666748, 113666822, 113676421,
                                                53963113),
                                     end   = c(113681827, 113681827, 113682211,
                                                54071569)),
                    strand = strand(c("+", "+", "+", "-")),
                    tx_id = c(29:31, 147L))
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

    gr <- GRanges(seqnames = "chr3",
                  ranges = IRanges(start=113677210, width=1),
                  strand = strand("+"))
    want <-
      GRanges(seqnames = factor("chr3", levels = seqlevels),
              ranges = IRanges(start=c(113677210, 113677210),
                               end  =c(113677385, 113682211)),
              strand = strand("+"),
              exon_id = 139:140)
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

    gr <- GRanges(seqnames = "chr3",
                  ranges = IRanges(start=113677210, width=1),
                  strand = strand("+"))
    want <-
      GRanges(seqnames = factor("chr3", levels = seqlevels),
              ranges = IRanges(start=c(113677210, 113677210),
                               end  =c(113677385, 113677477)),
              strand = strand("+"),
              cds_id = 116:117)
    seqinfo(want) <- seqinfo
    want <- GenomicFeatures:::.assignMetadataList(want, txdb)
    checkIdentical(cdsByOverlaps(txdb, gr), want)
}
