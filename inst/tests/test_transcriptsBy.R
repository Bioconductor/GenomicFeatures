test_transcriptsBy <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                     package="GenomicFeatures"))


    checkIdentical(transcriptsBy(txdb, "gene")[[1]],
                   GRanges(seqnames = "chr21_random",
                                  ranges   = IRanges(start=103280, end=164670),
                                  strand   = strand("-"),
                                  tx_name  = "uc002zka.1",
                                  tx_id    = 120L,
                                  gene_id  = "100132288"))
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
                           exon_id = c(1L,4L) ))
}



test_cdsBy <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                     package="GenomicFeatures"))

    checkIdentical(cdsBy(txdb, "gene")[[6]],
                   GRanges(seqnames = c("chr5","chr5"),
                           ranges = IRanges(start = c(258412,269844),
                                              end = c(259073,269964)),
                           strand = strand(c("-","-")),
                           cds_name = as.character(c(NA,NA)),
                           cds_id = c(53L,54L) ))
}
