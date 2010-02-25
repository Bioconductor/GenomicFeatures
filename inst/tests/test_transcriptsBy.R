test_transcriptsBy <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                     package="GenomicFeatures"))


    checkIdentical(transcriptsBy(txdb, "gene")[[1]],
                   GRanges(seqnames = "chr21_random",
                                  ranges = IRanges(start=103280, end=164670),
                                  strand = strand("-"),
                                  tx_name = "uc002zka.1"))
}



test_exonsBy <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                     package="GenomicFeatures"))

    checkIdentical(exonsBy(txdb, "tx")[[2]],
                   GRanges(seqnames = c("chr1","chr1"),
                           ranges = IRanges(start = c(1116,2476),
                                              end = c(2090,4272)),
                           strand = strand(c("+","+")),
                           exon_name = as.character(c(NA,NA)),
                           exon_rank = c(1L,2L),
                           tx_name = c("uc009vip.1","uc009vip.1"),
                           tx_chrom = c("chr1","chr1"),
                           tx_strand = c("+","+")
                           ))
}


test_cdsBy <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                     package="GenomicFeatures"))

    checkIdentical(cdsBy(txdb, "gene")[[6]],
                   GRanges(seqnames = c("chr5","chr5"),
                           ranges = IRanges(start = c(269844,258412),
                                              end = c(269964,259073)),
                           strand = strand(c("-","-")),
                           cds_name = as.character(c(NA,NA)),
                           exon_rank = c(2L,3L),
                           tx_name = c("uc003jam.1","uc003jam.1"),
                           tx_chrom = c("chr5","chr5"),
                           tx_strand = c("-","-")
                           ))
}
