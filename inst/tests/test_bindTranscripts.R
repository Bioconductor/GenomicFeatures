test_bindTranscripts <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                      package="GenomicFeatures"))

    checkException(bindTranscripts(txdb), silent = TRUE)
    checkException(bindTranscripts(txdb, IRanges()), silent = TRUE)

    ranges <- IRanges(start = c(1000, 1000, 20000, 30000),
                      end   = c(4000, 4000, 30000, 40000))
    chrom <- c("chr1", "chr1", "chr2", "chr2")
    strand <- strand(c("+", "-", "+", "-"))
    gf <- GenomicFeature(seqnames = chrom, ranges = ranges, strand = strand)
    want <- rd
    want[["tx_id"]] <- IntegerList(list(1:2, integer(), integer(), 4))
    want[["tx_name"]] <-
      CharacterList(list(c("uc001aaa.2","uc009vip.1"),
                         character(), character(), "uc002qvt.1"))
    checkIdentical(want, bindTranscripts(txdb, gf))
}

test_bindExons <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                     package="GenomicFeatures"))

    checkException(bindExons(txdb), silent = TRUE)
    checkException(bindExons(txdb, IRanges()), silent = TRUE)

    ranges <- IRanges(start = c(1000, 1000, 20000, 30000),
                      end   = c(4000, 4000, 30000, 40000))
    chrom <- c("chr1", "chr1", "chr2", "chr2")
    strand <- strand(c("+", "-", "+", "-"))
    gf <- GenomicFeature(seqnames = chrom, ranges = ranges, strand = strand)
    want <- rd
    want[["exon_id"]] <- IntegerList(list(c(1,2,4,3), integer(), integer(), 9:10))
    checkIdentical(want, bindExons(txdb, gf))
}

test_bindCDS <- function()
{
    suppressMessages(library(IRanges))
    suppressMessages(library(BSgenome))

    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                    package="GenomicFeatures"))

    checkException(bindCDS(txdb), silent = TRUE)
    checkException(bindCDS(txdb, IRanges()), silent = TRUE)

    ranges <- IRanges(start = c(1000, 1000, 20000, 30000),
                      end   = c(4000, 4000, 30000, 40000))
    chrom <- c("chr1", "chr1", "chr2", "chr2")
    strand <- strand(c("+", "-", "+", "-"))
    gf <- GenomicFeature(seqnames = chrom, ranges = ranges, strand = strand)
    want <- rd
    want[["cds_id"]] <- IntegerList(list(integer(), integer(), integer(), 1:2))
    checkIdentical(want, bindCDS(txdb, gf))
}
