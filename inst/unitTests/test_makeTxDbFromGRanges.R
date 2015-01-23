###

makeTxDbFromGRanges <- GenomicFeatures:::makeTxDbFromGRanges
format_txdb_dump <- GenomicFeatures:::.format_txdb_dump

test_makeTxDbFromGRanges_no_splicings <- function()
{
    ## 2 transcripts linked to a gene + 2 orphan transcripts
    gr_type   <- c( "mRNA", "mRNA",  "gene",  "mRNA", "mRNA")
    gr_ID     <- c( "tx1a",  "tx2", "gene1",  "tx1b",  "tx3")
    gr_Parent <- c("gene1",     NA,      NA, "gene1",     NA)
    gr_start  <- c(    151,    201,     101,     101,      5)
    gr_end    <- c(    160,    250,     160,     120,     35)
    gr <- GRanges("chr1", IRanges(gr_start, gr_end), strand="+")
    mcols(gr) <- DataFrame(type=gr_type, ID=gr_ID, Parent=gr_Parent)

    target_transcripts <- data.frame(
        tx_id=1:4,
        tx_name=c("tx3", "tx1b", "tx1a", "tx2"),
        tx_chrom="chr1",
        tx_strand="+",
        tx_start=c(5, 101, 151, 201),
        tx_end=c(35, 120, 160, 250)
    )
    target_genes <- data.frame(
        tx_id=2:3,
        gene_id="gene1"
    )
    target_chrominfo <- data.frame(
        chrom="chr1",
        length=NA,
        is_circular=NA
    )
    target_txdb_dump <- format_txdb_dump(transcripts=target_transcripts,
                                         splicings=NULL,
                                         genes=target_genes,
                                         chrominfo=target_chrominfo)

    current_txdb <- makeTxDbFromGRanges(gr)
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    mcols(gr)$Parent <- CharacterList(character(0))
    mcols(gr)$Parent[c(1L, 4L)] <- "gene1"
    current_txdb <- makeTxDbFromGRanges(gr)
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    ## 4 orphan transcripts
    gr4 <- gr[-3L]  # remove the gene

    target_txdb_dump <- format_txdb_dump(transcripts=target_transcripts,
                                         splicings=NULL,
                                         genes=NULL,
                                         chrominfo=target_chrominfo)

    current_txdb <- makeTxDbFromGRanges(gr4)
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    ## 2 orphan transcripts
    gr2 <- gr[c(2L, 5L)]

    target_transcripts2 <- target_transcripts[c(1,4), ]
    target_transcripts2$tx_id <- 1:2

    target_txdb_dump <- format_txdb_dump(transcripts=target_transcripts2,
                                         splicings=NULL,
                                         genes=NULL,
                                         chrominfo=target_chrominfo)

    current_txdb <- makeTxDbFromGRanges(gr2)
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    ## 0 transcript
    gr0 <- gr[0]

    target_txdb_dump <- format_txdb_dump(chrominfo=target_chrominfo)

    current_txdb <- makeTxDbFromGRanges(gr0)
    checkIdentical(target_txdb_dump, as.list(current_txdb))
}

test_makeTxDbFromGRanges_dup_tx_ids <- function()
{
}

