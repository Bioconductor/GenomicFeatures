test_transcriptsByRanges <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))

    checkException(transcriptsByRanges(txdb), silent = TRUE)
    checkException(transcriptsByRanges(txdb, IRanges()), silent = TRUE)
    checkException(transcriptsByRanges(txdb, RangedData(), columns = "bad"),
                   silent = TRUE)

    rd <- RangedData(ranges = IRanges(start=190000, end=280000),
                     space  = "chr5", strand = strand("-"))
    want <-
      RangedData(IRanges(start=257875, end=271297),
                 strand = strand("-"),
                 tx_id = 15L,
                 tx_name = "uc003jam.1",
                 space = "chr5")
    ranges(want) <- IRangesList(as.list(ranges(want)))
    checkIdentical(transcriptsByRanges(txdb, rd), want)
}

test_exonsByRanges <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))

    checkException(exonsByRanges(txdb), silent = TRUE)
    checkException(exonsByRanges(txdb, IRanges()), silent = TRUE)
    checkException(exonsByRanges(txdb, RangedData(), columns = "bad"),
                   silent = TRUE)

    rd <- RangedData(ranges = IRanges(start=190000, end=280000),
                     space  = "chr5", strand = strand("-"))
    want <-
      RangedData(IRanges(start=c(257875,269844,271208),
                         end  =c(259073,269974,271297)),
                 strand = strand(rep("-",3)),
                 exon_id = 77:79,
                 space = "chr5")
    ranges(want) <- IRangesList(as.list(ranges(want)))
    checkIdentical(exonsByRanges(txdb, rd), want)
}

test_cdsByRanges <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))

    checkException(cdsByRanges(txdb), silent = TRUE)
    checkException(cdsByRanges(txdb, IRanges()), silent = TRUE)
    checkException(cdsByRanges(txdb, RangedData(), columns = "bad"),
                   silent = TRUE)

    rd <- RangedData(ranges = IRanges(start=258412, end=269964),
                     space  = "chr5", strand = strand("-"))
    want <-
      RangedData(IRanges(start=c(258412,269844),
                         end  =c(259073,269964)),
                    strand = strand(rep("-",2)),
                    cds_id = 53:54,
                    space  = "chr5")
    ranges(want) <- IRangesList(as.list(ranges(want)))
    checkIdentical(cdsByRanges(txdb, rd), want)
}
