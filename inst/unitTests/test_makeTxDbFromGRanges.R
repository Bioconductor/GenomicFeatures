###

format_txdb_dump <- GenomicFeatures:::.format_txdb_dump

test_makeTxDbFromGRanges_empty <- function()
{
    ## The GFF3 way
    gr1 <- GRanges(type=character(0), ID=character(0))
    txdb1 <- makeTxDbFromGRanges(gr1)
    checkTrue(validObject(txdb1))
    checkTrue(all(elementNROWS(as.list(txdb1)) == 0L))
    checkIdentical(0L, length(seqinfo(txdb1)))
    checkIdentical(NA_character_, organism(txdb1))

    ## The GTF way
    gr2 <- GRanges(type=character(0),
                   gene_id=character(0), transcript_id=character(0))
    txdb2 <- makeTxDbFromGRanges(gr2)
    checkIdentical(as.list(txdb1), as.list(txdb2))
}

test_makeTxDbFromGRanges_no_exons_no_cds <- function()
{
    ## 2 transcripts linked to a gene + 2 orphan transcripts
    gene_ID    <- "gene1"
    gene_start <-     101
    gene_end   <-     160

    tx_ID      <- c( "tx1a", "tx2",  "tx1b", "tx3")
    tx_Parent  <- c("gene1",    NA, "gene1",    NA)
    tx_start   <- c(    151,   201,     101,     5)
    tx_end     <- c(    160,   250,     120,    35)

    gene_gr <- GRanges("chr1", IRanges(gene_start, gene_end), "+",
                       type="gene", ID=gene_ID, Parent=NA)
    tx_gr   <- GRanges("chr1", IRanges(tx_start, tx_end), "+",
                       type="mRNA", ID=tx_ID, Parent=tx_Parent)
    gr <- c(gene_gr, tx_gr)

    target_transcripts <- data.frame(
        tx_id=1:4,
        tx_name=c("tx3", "tx1b", "tx1a", "tx2"),
        tx_type="mRNA",
        tx_chrom="chr1",
        tx_strand="+",
        tx_start=c(5, 101, 151, 201),
        tx_end=c(35, 120, 160, 250)
    )
    target_splicings <- data.frame(
        tx_id=target_transcripts$tx_id,
        exon_rank=1,
        exon_id=1:4,
        exon_name=target_transcripts$tx_name,
        exon_chrom=target_transcripts$tx_chrom,
        exon_strand=target_transcripts$tx_strand,
        exon_start=target_transcripts$tx_start,
        exon_end=target_transcripts$tx_end
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
                                         splicings=target_splicings,
                                         genes=target_genes,
                                         chrominfo=target_chrominfo)

    current_txdb <- makeTxDbFromGRanges(gr)
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    ## make Parent a CharacterList instead of an atomic vector
    mcols(gr)$Parent <- CharacterList(character(0))
    mcols(gr)$Parent[c(2L, 4L)] <- "gene1"
    current_txdb <- makeTxDbFromGRanges(gr)
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    ## 4 orphan transcripts
    gr4 <- gr[mcols(gr)$type != "gene"]  # remove the gene

    target_txdb_dump <- format_txdb_dump(transcripts=target_transcripts,
                                         splicings=target_splicings,
                                         chrominfo=target_chrominfo)

    current_txdb <- makeTxDbFromGRanges(gr4)
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    ## 2 orphan transcripts
    gr2 <- gr[c(3L, 5L)]

    target_transcripts2 <- target_transcripts[c(1,4), ]
    target_transcripts2$tx_id <- 1:2
    target_splicings2 <- target_splicings[c(1,4), ]
    target_splicings2$tx_id <- 1:2
    target_splicings2$exon_id <- 1:2

    target_txdb_dump <- format_txdb_dump(transcripts=target_transcripts2,
                                         splicings=target_splicings2,
                                         chrominfo=target_chrominfo)

    current_txdb <- makeTxDbFromGRanges(gr2)
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    ## 1 gene with no children
    gr3 <- gr[-c(2L, 4L)]

    current_txdb <- makeTxDbFromGRanges(gr3)
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    ## 0 transcript
    gr0 <- gr[0]

    target_txdb_dump <- format_txdb_dump(chrominfo=target_chrominfo)

    current_txdb <- makeTxDbFromGRanges(gr0)
    checkIdentical(target_txdb_dump, as.list(current_txdb))   
}

test_makeTxDbFromGRanges_with_exons_and_cds <- function()
{
    ## 2 exons for tx1a, 1 exon for tx1b, 0 exon for tx2, 3 exons for tx3
    tx_ID       <- c( "tx1a", "tx1b", "tx2", "tx3")
    tx_start    <- c(    101,    145,   201,     5)
    tx_end      <- c(    160,    160,   250,    35)

    exon_Parent <- rep.int(tx_ID, c(2L, 1L, 0L, 3L))
    exon_start  <- c(101, 145, 145,  5, 18, 31)
    exon_end    <- c(134, 160, 160, 12, 28, 35)

    tx_gr   <- GRanges("chr1", IRanges(tx_start, tx_end), "+",
                       type="mRNA", ID=tx_ID, Parent=NA)
    exon_gr <- GRanges("chr1", IRanges(exon_start, exon_end), "+",
                       type="exon", ID=NA, Parent=exon_Parent)
    gr <- c(tx_gr, exon_gr)

    tx_oo <- c(4L, 1:3)
    target_transcripts <- data.frame(
        tx_id=1:4,
        tx_name=tx_ID[tx_oo],
        tx_type="mRNA",
        tx_chrom="chr1",
        tx_strand="+",
        tx_start=tx_start[tx_oo],
        tx_end=tx_end[tx_oo]
    )
    exon_oo <- c(4:6, 1:3)
    target_splicings <- data.frame(
        tx_id=c(rep(1:3, 3:1), 4),
        exon_rank=c(1:3, 1:2, 1, 1),
        exon_id=c(1:5, 5:6),
        exon_name=c(rep(NA, 6), "tx2"),
        exon_chrom="chr1",
        exon_strand="+",
        exon_start=c(exon_start[exon_oo], 201),
        exon_end=c(exon_end[exon_oo], 250)
    )
    target_chrominfo <- data.frame(
        chrom="chr1",
        length=NA,
        is_circular=NA
    )
    target_txdb_dump <- format_txdb_dump(transcripts=target_transcripts,
                                         splicings=target_splicings,
                                         chrominfo=target_chrominfo)

    current_txdb <- makeTxDbFromGRanges(gr)
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    ## naming exons thru ID
    exon_ID <- paste0("ex", seq_along(exon_gr))
    mcols(exon_gr)$ID <- exon_ID
    gr <- c(tx_gr, exon_gr)

    target_splicings$exon_id <- 1:7
    target_splicings$exon_name <- c(exon_ID[exon_oo], "tx2")
    target_txdb_dump <- format_txdb_dump(transcripts=target_transcripts,
                                         splicings=target_splicings,
                                         chrominfo=target_chrominfo)

    current_txdb <- makeTxDbFromGRanges(gr)
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    ## adding CDS
    cds_Parent  <- rep.int(tx_ID, c(2L, 1L, 0L, 2L))
    cds_start   <- c(122, 145, 149, 21, 31)
    cds_end     <- c(134, 152, 157, 28, 34)

    cds_gr  <- GRanges("chr1", IRanges(cds_start, cds_end), "+",
                       type="CDS", ID=NA, Parent=cds_Parent)
    gr <- c(gr, cds_gr)

    cds_oo <- c(NA, 4:5, 1:3)
    target_splicings <- cbind(
        target_splicings,
        cds_id=c(NA, 1:5, NA),
        cds_start=c(cds_start[cds_oo], NA),
        cds_end=c(cds_end[cds_oo], NA)
    )
    target_txdb_dump <- format_txdb_dump(transcripts=target_transcripts,
                                         splicings=target_splicings,
                                         chrominfo=target_chrominfo)

    current_txdb <- makeTxDbFromGRanges(gr)
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    ## with exons not separated by introns
    exon_end[c(1, 5)] <- c(144, 36)
    end(exon_gr) <- exon_end
    gr <- c(tx_gr, exon_gr)

    tx_oo <- 1:3
    target_transcripts <- data.frame(
        tx_id=1:3,
        tx_name=tx_ID[tx_oo],
        tx_type="mRNA",
        tx_chrom="chr1",
        tx_strand="+",
        tx_start=tx_start[tx_oo],
        tx_end=tx_end[tx_oo]
    )
    exon_oo <- 1:3
    target_splicings <- data.frame(
        tx_id=c(1, 1, 2, 3),
        exon_rank=c(1, 2, 1, 1),
        exon_id=1:4,
        exon_name=c(exon_ID[exon_oo], "tx2"),
        exon_chrom="chr1",
        exon_strand="+",
        exon_start=c(exon_start[exon_oo], 201),
        exon_end=c(exon_end[exon_oo], 250)
    )

    target_txdb_dump <- format_txdb_dump(transcripts=target_transcripts,
                                         splicings=target_splicings,
                                         chrominfo=target_chrominfo)

    current_txdb <- suppressWarnings(makeTxDbFromGRanges(gr))
    checkIdentical(target_txdb_dump, as.list(current_txdb))
}

test_makeTxDbFromGRanges_with_multiple_part_transcripts <- function()
{
    ## 1 gene with 2 transcripts + 1 gene with 1 transcript
    gene_ID     <- c("gene1", "gene2")
    gene_start  <- c(    101,       5)
    gene_end    <- c(    160,      35)

    tx_ID       <- c( "tx1a",  "tx1a",  "tx1b",   "tx2",   "tx2",   "tx2")
    tx_Parent   <- c("gene1", "gene1", "gene1", "gene2", "gene2", "gene2")
    tx_start    <- c(    145,     101,     145,      31,       5,      18)
    tx_end      <- c(    160,     134,     160,      35,      12,      28)

    gene_gr <- GRanges("chr1", IRanges(gene_start, gene_end), "+",
                       type="gene", ID=gene_ID, Parent=NA)
    tx_gr   <- GRanges("chr1", IRanges(tx_start, tx_end), "+",
                       type="mRNA", ID=tx_ID, Parent=tx_Parent)
    gr <- c(gene_gr, tx_gr)

    target_transcripts <- data.frame(
        tx_id=1:3,
        tx_name=c("tx2", "tx1a", "tx1b"),
        tx_type="mRNA",
        tx_chrom="chr1",
        tx_strand="+",
        tx_start=c(5, 101, 145),
        tx_end=c(35, 160, 160)
    )
    target_splicings <- data.frame(
        tx_id=target_transcripts$tx_id,
        exon_rank=1,
        exon_id=1:3,
        exon_name=target_transcripts$tx_name,
        exon_chrom=target_transcripts$tx_chrom,
        exon_strand=target_transcripts$tx_strand,
        exon_start=target_transcripts$tx_start,
        exon_end=target_transcripts$tx_end
    )
    target_genes <- data.frame(
        tx_id=1:3,
        gene_id=rep.int(c("gene2", "gene1"), c(1, 2))
    )
    target_chrominfo <- data.frame(
        chrom="chr1",
        length=NA,
        is_circular=NA
    )
    target_txdb_dump <- format_txdb_dump(transcripts=target_transcripts,
                                         splicings=target_splicings,
                                         genes=target_genes,
                                         chrominfo=target_chrominfo)

    current_txdb <- suppressWarnings(makeTxDbFromGRanges(gr))
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    ## adding 1 exon per transcript part
    exon_gr   <- GRanges("chr1", IRanges(tx_start, tx_end), "+",
                         type="exon", ID=NA, Parent=tx_ID)
    gr <- c(gr, exon_gr)

    exon_oo <- c(5:6, 4, 2, 1, 3)
    target_splicings <- data.frame(
        tx_id=rep.int(1:3, 3:1),
        exon_rank=c(1:3, 1:2, 1),
        exon_id=c(1:5, 5),
        exon_chrom="chr1",
        exon_strand="+",
        exon_start=tx_start[exon_oo],
        exon_end=tx_end[exon_oo]
    )
    target_txdb_dump <- format_txdb_dump(transcripts=target_transcripts,
                                         splicings=target_splicings,
                                         genes=target_genes,
                                         chrominfo=target_chrominfo)

    current_txdb <- suppressWarnings(makeTxDbFromGRanges(gr))
    checkIdentical(target_txdb_dump, as.list(current_txdb))
}

test_makeTxDbFromGRanges_with_exons_and_cds_as_gene_direct_children <-
    function()
{
    ## 1 gene with 2 transcripts + 1 gene with 3 exons as direct children +
    ## 1 gene with 2 CDS as direct children
    gene_ID     <- c("gene1", "gene2", "gene3")
    gene_start  <- c(    101,       5,     325)
    gene_end    <- c(    160,      35,     370)

    tx_ID       <- c( "tx1a",  "tx1b")
    tx_Parent   <- c("gene1", "gene1")
    tx_start    <- c(    101,     145)
    tx_end      <- c(    160,     160)

    exon_ID     <- c( "ex1",  "ex2",  "ex3",   "ex4",   "ex5",   "ex6")
    exon_Parent <- c("tx1a", "tx1a", "tx1b", "gene2", "gene2", "gene2")
    exon_start  <- c(   145,    101,    145,      31,       5,      18)
    exon_end    <- c(   160,    134,    160,      35,      12,      28)

    cds_Parent  <- c("gene3", "gene3")
    cds_start   <- c(    351,     325)
    cds_end     <- c(    370,     340)

    gene_gr <- GRanges("chr1", IRanges(gene_start, gene_end), "+",
                       type="gene", ID=gene_ID, Parent=NA)
    tx_gr   <- GRanges("chr1", IRanges(tx_start, tx_end), "+",
                       type="mRNA", ID=tx_ID, Parent=tx_Parent)
    exon_gr <- GRanges("chr1", IRanges(exon_start, exon_end), "+",
                       type="exon", ID=exon_ID, Parent=exon_Parent)
    cds_gr  <- GRanges("chr1", IRanges(cds_start, cds_end), "+",
                       type="CDS", ID=NA, Parent=cds_Parent)
    gr <- c(gene_gr, tx_gr, exon_gr, cds_gr)

    target_transcripts <- data.frame(
        tx_id=1:4,
        tx_name=c("gene2", "tx1a", "tx1b", "gene3"),
        tx_type=c("gene", "mRNA", "mRNA", "gene"),
        tx_chrom="chr1",
        tx_strand="+",
        tx_start=c(5, 101, 145, 325),
        tx_end=c(35, 160, 160, 370)
    )
    exon_oo <- c(5:6, 4, 2, 1, 3)
    cds_oo <- 2:1
    target_splicings <- data.frame(
        tx_id=rep.int(1:4, c(3:1, 2)),
        exon_rank=c(1:3, 1:2, 1, 1:2),
        exon_id=1:8,
        exon_name=c(exon_ID[exon_oo], NA, NA),
        exon_chrom="chr1",
        exon_strand="+",
        exon_start=c(exon_start[exon_oo], cds_start[cds_oo]),
        exon_end=c(exon_end[exon_oo], cds_end[cds_oo]),
        cds_id=c(rep.int(NA, length(exon_oo)), 1, 2),
        cds_start=c(rep.int(NA, length(exon_oo)), cds_start[cds_oo]),
        cds_end=c(rep.int(NA, length(exon_oo)), cds_end[cds_oo])
    )
    target_genes <- data.frame(
        tx_id=1:4,
        gene_id=rep.int(c("gene2", "gene1", "gene3"), c(1, 2, 1))
    )
    target_chrominfo <- data.frame(
        chrom="chr1",
        length=NA,
        is_circular=NA
    )
    target_txdb_dump <- format_txdb_dump(transcripts=target_transcripts,
                                         splicings=target_splicings,
                                         genes=target_genes,
                                         chrominfo=target_chrominfo)

    current_txdb <- makeTxDbFromGRanges(gr)
    checkIdentical(target_txdb_dump, as.list(current_txdb))

    ## Shuffling the input should not change anything.
    for (i in 1:5) {
        current_txdb <- makeTxDbFromGRanges(sample(gr))
        checkIdentical(target_txdb_dump, as.list(current_txdb))
    }
}

