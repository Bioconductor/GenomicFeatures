###

test_exonicParts <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                   package="GenomicFeatures"))

    exonic_parts <- exonicParts(txdb, linked.to.single.gene.only=TRUE)
    checkIdentical(length(exonic_parts), 653L)
    expected_mcolnames <- c("tx_id", "tx_name", "gene_id", "exon_id",
                            "exon_name", "exon_rank", "exonic_part")
    checkIdentical(expected_mcolnames, names(mcols(exonic_parts)))
    checkTrue(is.character(mcols(exonic_parts)$gene_id))
    checkTrue(is(mcols(exonic_parts)$tx_name, "CharacterList"))
    checkTrue(is.integer(mcols(exonic_parts)$exonic_part))

    exonic_parts <- exonicParts(txdb, linked.to.single.gene.only=FALSE)
    checkIdentical(length(exonic_parts), 660L)
    expected_mcolnames <- head(expected_mcolnames, n=-1L)
    checkIdentical(expected_mcolnames, names(mcols(exonic_parts)))
}

