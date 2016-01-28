###

test_makeTxDbFromBiomart <- function()
{
    ## We test makeTxDbFromBiomart() on the following transcripts:
    TARGET_TX_NAME <- c(
        ## Coding transcripts.
        "ENST00000013894",  #  6 exons on + strand, CDS'es on exons 1:3
        "ENST00000268655",  #  3 exons on + strand, CDS'es on all exons
        "ENST00000313243",  # 14 exons on - strand, CDS'es on exons 2:14
        "ENST00000435657",  # 30 exons on - strand, CDS'es on exons 2:30
        ## Non-coding transcripts.
        "ENST00000384428",  #  1 exon  on - strand, no CDS
        "ENST00000478783"   #  2 exons on + strand, no CDS
    )
    TARGET_NEXONS_PER_TX <- c(6L, 3L, 14L, 30L, 1L, 2L)
    TARGET_TX_STRAND <- c("+", "+", "-", "-", "-", "+")
    TARGET_CDS2EXON <- list(1:3, 1:3, 2:14, 2:30)
    TARGET_GENE <- c("ENSG00000011198", "ENSG00000103343",
                     "ENSG00000111837", "ENSG00000231116",
                     "ENSG00000207157", "ENSG00000067646")

    current_txdb <- makeTxDbFromBiomart(transcript_ids=TARGET_TX_NAME,
                                        circ_seqs="MT")
    checkTrue(validObject(current_txdb))

    ## Extract transcripts and re-order them as in TARGET_TX_NAME.
    current_tx <- transcripts(current_txdb, columns=c("tx_name", "gene_id"))
    checkIdentical(length(TARGET_TX_NAME), length(current_tx))
    current_tx_name <- mcols(current_tx)$tx_name
    target2current <- match(TARGET_TX_NAME, current_tx_name)
    current_tx <- current_tx[target2current]

    ## Check transcript name.
    current_tx_name <- mcols(current_tx)$tx_name
    checkIdentical(TARGET_TX_NAME, current_tx_name)

    ## Check strand.
    current_tx_strand <- as.character(strand(current_tx))
    checkIdentical(TARGET_TX_STRAND, current_tx_strand)

    ## Check gene.
    current_gene <- as.character(mcols(current_tx)$gene_id)
    checkIdentical(TARGET_GENE, current_gene)
 
    ## Check nb of exons per transcript.
    ex_by_tx <- exonsBy(current_txdb, by="tx", use.names=TRUE)
    checkTrue(setequal(TARGET_TX_NAME, names(ex_by_tx)))
    nexons_per_tx <- elementNROWS(ex_by_tx)[TARGET_TX_NAME]
    checkIdentical(TARGET_NEXONS_PER_TX, as.integer(nexons_per_tx))

    ## Check CDS'es.
    tx_names_with_cds <- head(TARGET_TX_NAME, n=-2)
    cds_by_tx <- cdsBy(current_txdb, by="tx", use.names=TRUE)
    checkTrue(setequal(tx_names_with_cds, names(cds_by_tx)))
    cds_by_tx <- cds_by_tx[tx_names_with_cds]
    current_cds2exon <- mcols(unlist(cds_by_tx, use.names=FALSE))$exon_rank
    checkIdentical(unlist(TARGET_CDS2EXON), current_cds2exon)
  
    ## Compare with a precomputed TxDb object.
    target_file <- system.file("extdata",
                               "Biomart_Ensembl_sample.sqlite",
                               package="GenomicFeatures")
    target_txdb <- loadDb(target_file)
    checkTrue(GenomicFeatures:::compareTxDbs(target_txdb, current_txdb))
}

