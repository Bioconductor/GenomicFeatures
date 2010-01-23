## TODO: find out why these are not the SAME!
## test_transcripts <- function()
## {    
##     txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
##                                       package="GenomicFeatures"))
##     suppressMessages(library(IRanges))
##     vals <- list(tx_chrom = c("chr1", "chr5"), tx_strand = "-")
    
##     wantRanges = IRanges(start = c(4269,257875),
##                          end   = c(6628,271297))
##     want<-RangedData(ranges  = wantRanges,
##                      strand  = factor(rep("-",2),
##                                      levels=c("-","+","*")),
##                      exon_id = IntegerList("69"=c(79,78,77),"110"=c(8,7,6,5)),
##                      space   = c("chr1","chr5"),
##                      tx_id   = c(110L,69L),
##                      tx_name = c("uc009vis.1","uc003jam.1"))

##     checkEquals(want, transcripts(txdb, vals, columns = c("tx_id","tx_name","exon_id")) )
## }

## transcripts(txdb, vals, columns = c("tx_id","tx_name","exon_id"))

test_exons <- function()
{
    txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))
    suppressMessages(library(IRanges))
    vals <- list(tx_chrom = c("chr1", "chr5"), tx_strand = "-")
    
    wantRanges = IRanges(start = c(4269,4833,5659,6470,257875,269844,271208),
                         end   = c(4692,4901,5805,6628,259073,269974,271297))

    
    want<-RangedData(ranges  = wantRanges,
                     strand  = factor(rep("-",7),
                                     levels=c("-","+","*")),
                     space   = c(rep("chr1",4),rep("chr5",3)),
                     exon_id   = c(5,6,7,8,77,78,79))

    checkEquals(want, exons(txdb, vals))    
}

