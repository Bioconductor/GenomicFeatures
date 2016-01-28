###

test_disjointExons <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                   package="GenomicFeatures"))
    de <- disjointExons(txdb, FALSE, TRUE)
    checkTrue(is(de$gene_id, "CharacterList"))
    checkTrue(is(de$tx_name, "CharacterList"))
    checkTrue(is(de$exonic_part, "integer"))
    checkIdentical(sum(elementNROWS(de$gene_id) == 1), 653L)
    checkIdentical(sum(elementNROWS(de$gene_id) == 2), 0L)

    de <- disjointExons(txdb, TRUE, FALSE)
    checkTrue(all(names(mcols(de)) %in% c("gene_id", "exonic_part")))
    checkIdentical(sum(elementNROWS(de$gene_id) == 1), 653L)
    checkIdentical(sum(elementNROWS(de$gene_id) == 2), 0L)
}

