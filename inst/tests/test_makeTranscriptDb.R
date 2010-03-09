###

test_makeTranscriptDb <- function()
{
    ## AN UNREALISTIC EDGE CASE
    ## ------------------------
    transcripts0 <- data.frame(
                        tx_id=c(26L, 5L, 11L),
                        tx_chrom=c("chr1", "chr2", "chr2"),
                        tx_strand=c("+", "-", "-"),
                        tx_start=c(1L, 16844685L, 16844685L),
                        tx_end=c(100L, 16844760L, 16844760L))
    splicings0 <-  data.frame(
                        tx_id=c(26L, 5L, 26L, 11L),
                        exon_rank=c(2L, 1L, 1L, 1L),
                        exon_start=c(1L, 16844685L, 1L, 16844685L),
                        exon_end=c(100L, 16844760L, 100L, 16844760L))
    txdb0 <- makeTranscriptDb(transcripts0, splicings0)

    transcripts1 <- data.frame(
                        tx_id=transcripts0$tx_id,
                        tx_name=as.character(c(NA, NA, NA)),
                        tx_chrom=transcripts0$tx_chrom,
                        tx_strand=transcripts0$tx_strand,
                        tx_start=transcripts0$tx_start,
                        tx_end=transcripts0$tx_end,
                        stringsAsFactors=FALSE)
    checkIdentical(as.list(txdb0)$transcripts, transcripts1)

    splicings1 <- splicings0[c(3L, 1L, 2L, 4L), ]
    splicings1_exons <- data.frame(
                        exon_id=c(1L, 1L, 2L, 2L),
                        exon_name=as.character(c(NA, NA, NA, NA)),
                        exon_chrom=factor(c("chr1", "chr1", "chr2", "chr2")),
                        exon_strand=factor(c("+", "+", "-", "-")),
                        stringsAsFactors=FALSE)
    splicings1_cds <- data.frame(
                        cds_id=NA_integer_,
                        cds_start=NA_integer_,
                        cds_end=NA_integer_)[rep.int(1L, 4L), ]
    splicings1 <- cbind(splicings1[1:2],
                        splicings1_exons,
                        splicings1[3:4],
                        splicings1_cds)
    row.names(splicings1) <- NULL
    checkIdentical(as.list(txdb0)$splicings, splicings1)

    ## WITH REAL DATA
    ## -------------- 
    ## want
    txdb0_file <- system.file(
                      "extdata",
                      "Biomart_Ensembl_sample.sqlite",
                      package="GenomicFeatures")
    txdb0 <- loadFeatures(txdb0_file)

    ## get
    txdb1 <- do.call(makeTranscriptDb, as.list(txdb0))

    ## compare
    ok <- GenomicFeatures:::compareTranscriptDbs(txdb1, txdb0)
    checkTrue(ok)
}

test_makeTranscriptDbFromUCSCTxTable <- function()
{
    ## want
    txdb0_file <- system.file(
                      "extdata",
                      "UCSC_knownGene_sample.sqlite",
                      package="GenomicFeatures"
                  )
    txdb0 <- loadFeatures(txdb0_file)

    ## get
    knownGene_sample_file <- system.file(
                                 "extdata",
                                 "UCSC_knownGene_sample.rda",
                                 package="GenomicFeatures"
                             )
    load(knownGene_sample_file)
    knownToLocusLink_sample_file <- system.file(
                                        "extdata",
                                        "UCSC_knownToLocusLink_sample.rda",
                                        package="GenomicFeatures"
                                    )
    load(knownToLocusLink_sample_file)
    genes <- data.frame(tx_name=UCSC_knownToLocusLink_sample$name,
                        gene_id=UCSC_knownToLocusLink_sample$value)
    txdb1 <- GenomicFeatures:::.makeTranscriptDbFromUCSCTxTable(
                 UCSC_knownGene_sample, genes,
                 "hg18", "knownGene", "Entrez Gene ID", FALSE)

    ## compare
    ok <- GenomicFeatures:::compareTranscriptDbs(txdb1, txdb0)
    checkTrue(ok)
}

test_makeTranscriptDbFromBiomart <- function()
{
    ## want
    txdb0_file <- system.file(
                      "extdata",
                      "Biomart_Ensembl_sample.sqlite",
                      package="GenomicFeatures")
    txdb0 <- loadFeatures(txdb0_file)

    ## get
    transcript_ids <- c(
         "ENST00000400839",
         "ENST00000400840",
         "ENST00000478783",
         "ENST00000435657",
         "ENST00000268655",
         "ENST00000313243",
         "ENST00000341724"
    )
    txdb1 <- makeTranscriptDbFromBiomart(transcript_ids=transcript_ids)

    ## compare
    ok <- GenomicFeatures:::compareTranscriptDbs(txdb1, txdb0)
    checkTrue(ok)
}

