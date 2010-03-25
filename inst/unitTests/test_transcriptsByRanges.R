test_transcriptsByRanges <- function()
{
    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))
    seqlengths <- seqlengths(txdb)
    seqlevels <- names(seqlengths)

    checkException(transcriptsByRanges(txdb), silent = TRUE)
    checkException(transcriptsByRanges(txdb, IRanges()), silent = TRUE)
    checkException(transcriptsByRanges(txdb, GRanges(), columns = "bad"),
                   silent = TRUE)

    gr <- GRanges(seqnames = "chr5",
                  ranges = IRanges(start=190000, end=280000),
                  strand = strand("-"))
    want <-
      GRanges(seqnames = factor("chr5", levels = seqlevels),
              ranges = IRanges(start=257875, end=271297),
              strand = strand("-"),
              tx_id = 15L,
              tx_name = "uc003jam.1",
              seqlengths = seqlengths)
    checkIdentical(transcriptsByRanges(txdb, gr), want)

    ranges <- IRanges(start = c(1000, 1000, 20000, 30000),
                      end   = c(4000, 4000, 30000, 40000))
    chrom <- c("chr1", "chr1", "chr2", "chr2")
    strand <- strand(c("+", "-", "+", "-"))
    gr <- GRanges(seqnames = chrom, ranges = ranges, strand = strand)
    want <- GRanges(seqnames =
                    factor(c("chr1", "chr1", "chr2"), levels = seqlevels),
                    ranges = IRanges(start = c(1116, 1116, 31608),
                                     end   = c(4121, 4272, 36385)),
                    strand = strand(c("+", "+", "-")),
                    tx_id = c(1L,2L,4L),
                    seqlengths = seqlengths)
    checkIdentical(want, transcriptsByRanges(txdb, gr, columns="tx_id"))
}

test_exonsByRanges <- function()
{
    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))
    seqlengths <- seqlengths(txdb)
    seqlevels <- names(seqlengths)

    checkException(exonsByRanges(txdb), silent = TRUE)
    checkException(exonsByRanges(txdb, IRanges()), silent = TRUE)
    checkException(exonsByRanges(txdb, GRanges(), columns = "bad"),
                   silent = TRUE)

    gr <- GRanges(seqnames = "chr5",
                  ranges = IRanges(start=190000, end=280000),
                  strand = strand("-"))
    want <-
      GRanges(seqnames = factor(rep("chr5",3), levels = seqlevels),
              ranges = IRanges(start=c(257875,269844,271208),
                               end  =c(259073,269974,271297)),
              strand = strand(rep("-",3)),
              exon_id = 77:79,
              seqlengths = seqlengths)
    checkIdentical(exonsByRanges(txdb, gr), want)
}

test_cdsByRanges <- function()
{
    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))
    seqlengths <- seqlengths(txdb)
    seqlevels <- names(seqlengths)

    checkException(cdsByRanges(txdb), silent = TRUE)
    checkException(cdsByRanges(txdb, IRanges()), silent = TRUE)
    checkException(cdsByRanges(txdb, GRanges(), columns = "bad"),
                   silent = TRUE)

    gr <- GRanges(seqnames  = "chr5",
                  ranges = IRanges(start=258412, end=269964),
                  strand = strand("-"))
    want <-
      GRanges(seqnames = factor(rep("chr5",2), levels = seqlevels),
              ranges = IRanges(start=c(258412,269844),
                               end  =c(259073,269964)),
              strand = strand(rep("-",2)),
              cds_id = 53:54,
              seqlengths = seqlengths)
    checkIdentical(cdsByRanges(txdb, gr), want)
}
