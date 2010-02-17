test_transcriptsByRanges <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))

    checkException(transcriptsByRanges(txdb), silent = TRUE)
    checkException(transcriptsByRanges(txdb, IRanges()), silent = TRUE)
    checkException(transcriptsByRanges(txdb, GenomicFeature(), columns = "bad"),
                   silent = TRUE)

    gf <- GenomicFeature(seqnames = "chr5",
                         ranges = IRanges(start=190000, end=280000),
                         strand = strand("-"))
    want <-
      split(GenomicFeature(seqnames = "chr5",
                           ranges = IRanges(start=257875, end=271297),
                           strand = strand("-"),
                           tx_id = 15L,
                           tx_name = "uc003jam.1"))
    ranges(want) <- IRangesList(as.list(ranges(want)))
    checkIdentical(transcriptsByRanges(txdb, gf), want)

    ranges <- IRanges(start = c(1000, 1000, 20000, 30000),
                      end   = c(4000, 4000, 30000, 40000))
    chrom <- c("chr1", "chr1", "chr2", "chr2")
    strand <- strand(c("+", "-", "+", "-"))
    gf <- GenomicFeature(seqnames = chrom, ranges = ranges, strand = strand)
    want <- GenomicFeature(space = c("chr1", "chr1", "chr2"),
                           ranges = IRanges(start = c(1116, 1116, 31608),
                                            end   = c(4121, 4272, 36385)),
                           strand = strand(c("+", "+", "-")),
                           tx_id = c(1L,2L,4L))
    checkIdentical(want, transcriptsByRanges(txdb, gf, columns="tx_id"))
}

test_exonsByRanges <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))

    checkException(exonsByRanges(txdb), silent = TRUE)
    checkException(exonsByRanges(txdb, IRanges()), silent = TRUE)
    checkException(exonsByRanges(txdb, GenomicFeature(), columns = "bad"),
                   silent = TRUE)

    gf <- GenomicFeature(seqnames = "chr5",
                         ranges = IRanges(start=190000, end=280000),
                         strand = strand("-"))
    want <-
      GenomicFeature(seqnames = "chr5",
                     ranges = IRanges(start=c(257875,269844,271208),
                                      end  =c(259073,269974,271297)),
                     strand = strand(rep("-",3)),
                     exon_id = 77:79)
    ranges(want) <- IRangesList(as.list(ranges(want)))
    checkIdentical(exonsByRanges(txdb, gf), want)
}

test_cdsByRanges <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))

    checkException(cdsByRanges(txdb), silent = TRUE)
    checkException(cdsByRanges(txdb, IRanges()), silent = TRUE)
    checkException(cdsByRanges(txdb, GenomicFeature(), columns = "bad"),
                   silent = TRUE)

    gf <- GenomicFeature(seqnames  = "chr5",
                         ranges = IRanges(start=258412, end=269964),
                         strand = strand("-"))
    want <-
      GenomicFeature(seqnames = "chr5",
                     ranges = IRanges(start=c(258412,269844),
                                      end  =c(259073,269964)),
                     strand = strand(rep("-",2)),
                     cds_id = 53:54))
    ranges(want) <- IRangesList(as.list(ranges(want)))
    checkIdentical(cdsByRanges(txdb, gf), want)
}
