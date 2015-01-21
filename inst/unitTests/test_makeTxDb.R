###

test_makeTxDb <- function()
{
    ## A TOY CASE
    ## ----------
    transcripts0 <- data.frame(
                        tx_id=c(26L, 5L, 11L),
                        tx_name=c("A", "B", "C"),
                        tx_chrom=c("chr1", "chr2", "chr2"),
                        tx_strand=c("+", "-", "-"),
                        tx_start=c(1L, 16844685L, 16844685L),
                        tx_end=c(100L, 16844760L, 16844760L))
    splicings0 <-  data.frame(
                        tx_id=c(26L, 5L, 26L, 11L),
                        exon_rank=c(2L, 1L, 1L, 1L),
                        exon_start=c(1L, 16844685L, 1L, 16844685L),
                        exon_end=c(100L, 16844760L, 100L, 16844760L))

    suppressWarnings(txdb0a <- makeTxDb(transcripts0, splicings0,
                                                reassign.ids=TRUE))
    suppressWarnings(txdb0b <- makeTxDb(transcripts0, splicings0))

    ## Check the transcripts.
    transcripts1 <- data.frame(
                        tx_id=transcripts0$tx_id,
                        tx_name=as.character(transcripts0$tx_name),
                        tx_chrom=transcripts0$tx_chrom,
                        tx_strand=transcripts0$tx_strand,
                        tx_start=transcripts0$tx_start,
                        tx_end=transcripts0$tx_end,
                        stringsAsFactors=FALSE)
    oo <- order(transcripts1$tx_id)
    transcripts1b <- transcripts1[oo, ]
    rownames(transcripts1b) <- NULL
    transcripts1a <- transcripts1
    transcripts1a$tx_id <- 1:3 
    checkIdentical(as.list(txdb0a)$transcripts, transcripts1a)
    checkIdentical(as.list(txdb0b)$transcripts, transcripts1b)

    ## Check the splicings.
    rowmap <- match(splicings0$tx_id, transcripts0$tx_id)
    splicings1 <- data.frame(
                        tx_id=splicings0$tx_id,
                        exon_rank=splicings0$exon_rank,
                        exon_id=c(1L, 2L, 1L, 2L),
                        exon_chrom=transcripts0$tx_chrom[rowmap],
                        exon_strand=transcripts0$tx_strand[rowmap],
                        exon_start=splicings0$exon_start,
                        exon_end=splicings0$exon_end,
                        cds_id=as.integer(c(NA, NA, NA, NA)),
                        cds_start=as.integer(c(NA, NA, NA, NA)),
                        cds_end=as.integer(c(NA, NA, NA, NA)),
                        stringsAsFactors=FALSE)
    splicings1a <- splicings1
    splicings1a$tx_id <- c(1L, 2L, 1L, 3L)
    oo <- order(splicings1a$tx_id, splicings1a$exon_rank)
    splicings1a <- splicings1a[oo, ]
    rownames(splicings1a) <- NULL
    splicings1b <- splicings1
    oo <- order(splicings1b$tx_id, splicings1b$exon_rank)
    splicings1b <- splicings1b[oo, ]
    rownames(splicings1b) <- NULL
    checkIdentical(as.list(txdb0a)$splicings, splicings1a)
    checkIdentical(as.list(txdb0b)$splicings, splicings1b)

    ## WITH REAL DATA
    ## -------------- 
    ## want
    txdb0_file <- system.file(
                      "extdata",
                      "Biomart_Ensembl_sample.sqlite",
                      package="GenomicFeatures")
    txdb0 <- loadDb(txdb0_file)

    ## get
    txdb1 <- do.call(makeTxDb, as.list(txdb0))

    ## compare
    ok <- GenomicFeatures:::compareTxDbs(txdb1, txdb0)
    checkTrue(ok)
}

