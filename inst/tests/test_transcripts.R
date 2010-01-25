test_transcripts <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                     package="GenomicFeatures"))

    checkException(transcripts(data.frame()), silent = TRUE)
    checkException(transcripts(txdb, vals = list("bad" = 1:10)), silent = TRUE)
    checkException(transcripts(txdb, columns = "bad"), silent = TRUE)

    checkIdentical(transcripts(txdb, list("tx_id" = 3)),
                   RangedData(IRanges(start=4269, end=6628),
                              strand = strand("-"),
                              tx_id = 3L,
                              tx_name = "uc009vis.1",
                              space = "chr1"))

    vals <- list(tx_chrom = c("chr1", "chr5"), tx_strand = "-")
    wantRanges <- IRanges(start = c(4269, 257875),
                          end   = c(6628, 271297))
    want <- RangedData(ranges  = wantRanges,
                       strand  = strand(rep("-", 2)),
                       tx_id   = c(3L, 15L),
                       tx_name = c("uc009vis.1", "uc003jam.1"),
                       exon_id = IntegerList("3"=c(8,7,6,5), "15"=c(79,78,77)),
                       space   = c("chr1", "chr5"))
    want[["exon_id"]] <- IntegerList("3"=c(8,7,6,5), "15"=c(79,78,77))
    
    checkIdentical(want,
                   transcripts(txdb, vals, columns = c("tx_id","tx_name","exon_id")))
}

test_exons <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))

    checkException(exons(data.frame()), silent = TRUE)
    checkException(exons(txdb, vals = list("bad" = 1:10)), silent = TRUE)

    checkIdentical(exons(txdb, list("exon_id" = 1)),
                   RangedData(IRanges(start=1116, end=2090),
                              strand = strand("+"),
                              exon_id = 1L,
                              space = "chr1"))

    vals <- list(exon_chrom = c("chr1", "chr5"), exon_strand = "-")
    wantRanges <- IRanges(start = c(4269,4833,5659,6470,257875,269844,271208),
                          end   = c(4692,4901,5805,6628,259073,269974,271297))
    want <- RangedData(ranges  = wantRanges,
                       strand  = strand(rep("-", 7)),
                       exon_id   = c(5L, 6L, 7L, 8L, 77L, 78L, 79L),
                       space   = c(rep("chr1", 4), rep("chr5", 3)))

    checkIdentical(want, exons(txdb, vals))    
}
