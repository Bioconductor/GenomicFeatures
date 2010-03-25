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

