###

test_transcriptsByOverlaps <- function()
{
    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))

    checkException(transcriptsByOverlaps(txdb), silent = TRUE)
    checkException(transcriptsByOverlaps(txdb, IRanges()), silent = TRUE)
    checkException(transcriptsByOverlaps(txdb, GRanges(), columns = "bad"),
                   silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    gr <- GRanges(seqnames = "chr5",
                  ranges = IRanges(start=190000, end=280000),
                  strand = strand("-"))
    want <-
      GRanges(seqnames = factor("chr5", levels = seqlevels),
              ranges = IRanges(start=257875, end=271297),
              strand = strand("-"),
              tx_id = 15L,
              tx_name = "uc003jam.1")
    seqinfo(want) <- seqinfo
    metadata(want)[[1]] <- DataFrame(metadata(txdb))
    checkIdentical(transcriptsByOverlaps(txdb, gr), want)

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
                    tx_id = c(1L,2L,4L))
    seqinfo(want) <- seqinfo
    metadata(want)[[1]] <- DataFrame(metadata(txdb))
    checkIdentical(transcriptsByOverlaps(txdb, gr, columns="tx_id"), want)
}

test_exonsByOverlaps <- function()
{
    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))

    checkException(exonsByOverlaps(txdb), silent = TRUE)
    checkException(exonsByOverlaps(txdb, IRanges()), silent = TRUE)
    checkException(exonsByOverlaps(txdb, GRanges(), columns = "bad"),
                   silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    gr <- GRanges(seqnames = "chr5",
                  ranges = IRanges(start=190000, end=280000),
                  strand = strand("-"))
    want <-
      GRanges(seqnames = factor(rep("chr5",3), levels = seqlevels),
              ranges = IRanges(start=c(257875,269844,271208),
                               end  =c(259073,269974,271297)),
              strand = strand(rep("-",3)),
              exon_id = 77:79)
    seqinfo(want) <- seqinfo
    metadata(want)[[1]] <- DataFrame(metadata(txdb))
    checkIdentical(exonsByOverlaps(txdb, gr), want)
}

test_cdsByOverlaps <- function()
{
    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))

    checkException(cdsByOverlaps(txdb), silent = TRUE)
    checkException(cdsByOverlaps(txdb, IRanges()), silent = TRUE)
    checkException(cdsByOverlaps(txdb, GRanges(), columns = "bad"),
                   silent = TRUE)

    seqinfo <- seqinfo(txdb)
    seqlevels <- seqlevels(seqinfo)

    gr <- GRanges(seqnames  = "chr5",
                  ranges = IRanges(start=258412, end=269964),
                  strand = strand("-"))
    want <-
      GRanges(seqnames = factor(rep("chr5",2), levels = seqlevels),
              ranges = IRanges(start=c(258412,269844),
                               end  =c(259073,269964)),
              strand = strand(rep("-",2)),
              cds_id = 53:54)
    seqinfo(want) <- seqinfo
    metadata(want)[[1]] <- DataFrame(metadata(txdb))
    checkIdentical(cdsByOverlaps(txdb, gr), want)
}
