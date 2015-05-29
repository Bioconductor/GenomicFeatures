library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
cdsbytx <- cdsBy(txdb, "tx", use.names=TRUE)[1:3]

test_mapToTranscripts <- function()
{
    ## empty
    ans1 <- mapToTranscripts(GRanges(), GRanges())
    ans2 <- mapToTranscripts(GRanges(), GRangesList())
    checkTrue(length(ans1) == 0L)
    checkTrue(length(ans2) == 0L)
    checkIdentical(names(mcols(ans1)), c("xHits", "transcriptsHits"))
    checkIdentical(names(mcols(ans2)), c("xHits", "transcriptsHits"))

    ## strand
    x <- GRanges("chr2L", IRanges(c(7500, 8400, 9000), 
                 width=200, names=LETTERS[1:3]), strand="+")
    align <- cdsbytx

    ans <- mapToTranscripts(x, cdsbytx, FALSE, FALSE)
    checkTrue(length(ans) == 3L)
    checkIdentical(names(ans), c("B", "B", "C"))
    checkTrue(all(width(ans) == 200L))
    checkIdentical(mcols(ans)$xHits, c(2L, 2L, 3L))
    checkIdentical(mcols(ans)$transcriptsHits, c(1L, 3L, 2L))

    strand(align[[2]][1]) <- "-"
    checkException(mapToTranscripts(x, align, FALSE, FALSE), silent=TRUE) 

    strand(align[[2]]) <- "-"
    ans <- mapToTranscripts(x, align, ignore.strand=FALSE) 
    checkIdentical(mcols(ans)$xHits, c(2L, 2L))
    checkIdentical(mcols(ans)$transcriptsHits, c(1L, 3L))
    checkIdentical(names(ans), c("B", "B"))
    checkIdentical(seqlevels(ans), c("FBtr0300689", "FBtr0330654"))

    ## TxDb
    x <- GRanges("chr2L", IRanges(c(7500, 8400, 9000), 
                 width=200, names=LETTERS[1:3]))
    checkException(mapToTranscripts(x, txdb, extractor.fun="foo"), silent=TRUE)
    ans <- mapToTranscripts(x, txdb, extractor.fun=exonsBy, by="tx")
    checkIdentical(names(ans), c("B", "B", rep("C", 3)))
    checkIdentical(as.character(seqnames(ans)), as.character(c(1, 3, 1, 2, 3))) 
    ans <- mapToTranscripts(x, txdb, extractor.fun=exonsBy, by="gene")
    checkIdentical(mcols(ans)$transcriptsHits, rep(4179L, 5))
    checkIdentical(seqlevels(ans), "FBgn0031208")
}

test_mapFromTranscripts <- function()
{
    ## strand
    x <- GRanges("chr1", IRanges(rep(10, 5), width=2), strand="+")
    names(x) <- c("A", "B", "C", "foo", "A")
    gr <- GRanges("chr1", IRanges(c(10, 100), width=50), strand="+")
    align = GRangesList(gr, gr, gr, gr)
    names(align) <- c("C", "bar", "C", "A") 

    ans <- mapFromTranscripts(x, align)
    checkTrue(length(ans) == 4L)
    checkIdentical(names(ans), c("A", "C", "C", "A"))
    checkTrue(all(width(ans) == 2L))
    checkIdentical(mcols(ans)$xHits, c(1L, 3L, 3L, 5L))
    checkIdentical(mcols(ans)$transcriptsHits, c(4L, 1L, 3L, 4L))
    ans <- mapFromTranscripts(x, unlist(align))
    checkTrue(length(ans) == 8L)
    checkIdentical(names(ans), c("A", "A", rep("C", 4), "A", "A"))
    expected <- c(7L, 8L, 1L, 2L, 5L, 6L, 7L, 8L)
    checkIdentical(mcols(ans)$transcriptsHits, expected) 

    strand(align[[1]][1]) <- "-"
    checkException(mapFromTranscripts(x, align, FALSE), silent=TRUE) 

    strand(align[[1]]) <- "-"
    ans <- mapFromTranscripts(x, align, FALSE) 
    checkIdentical(mcols(ans)$transcriptsHits, c(4L, 3L, 4L))
    ans <- mapFromTranscripts(x, align, ignore.strand=TRUE) 
    checkIdentical(mcols(ans)$transcriptsHits, c(4L, 1L, 3L, 4L))
    checkIdentical(seqlevels(ans), "chr1")
}

test_mapToTranscripts_range_order <- function()
{
    x <- GRanges("chrA", IRanges(43522349, width=1), strand="+")
    align1 <- GRangesList(GRanges("chrA", 
        IRanges(c(43522244, 43528406),
                c(43524145, 43528644)), strand="+"))
    align2 <- GRangesList(GRanges("chrA", 
        IRanges(c(43528406, 43522244),
                c(43528644, 43524145)), strand="+"))
    names(align1) <- names(align2) <- LETTERS[seq_along(align1)]

    ## "+" strand
    ## smallest range first 
    ans <- mapToTranscripts(x, align1, FALSE)
    checkTrue(start(ans) == 106L)
    ## largest range first 
    ans <- mapToTranscripts(x, align2, FALSE)
    checkTrue(start(ans) == 106L)

    ## "-" strand
    strand(x) <- "-"
    strand(align1) <- "-"
    strand(align2) <- "-"
    ## smallest range first
    ans <- mapToTranscripts(x, align1, FALSE)
    checkTrue(start(ans) == 2036L)
    ## largest range first
    ans <- mapToTranscripts(x, align2, FALSE)
    checkTrue(start(ans) == 2036L)
}

test_pmapToTranscripts <- function()
{
    ## empty
    ans1 <- pmapToTranscripts(GRanges(), GRanges())
    ans2 <- pmapToTranscripts(GRanges(), GRangesList())
    checkTrue(length(ans1) == 0L)
    checkTrue(length(ans2) == 0L)
    checkIdentical(names(mcols(ans1)), character(0))
    checkIdentical(names(mcols(ans2)), character(0))

    ## length
    x <- GRanges("chr1", IRanges(1, width=1))
    y <- GRanges("chr1", IRanges(1:5, width=1))
    checkException(pmapToTranscripts(x, y), silent=TRUE)
    checkException(pmapToTranscripts(x, GRangesList(y, y)), silent=TRUE)

    ## strand
    x <- GRanges("chr1", IRanges(c(6, 16, 1), width=1, names=LETTERS[1:3]))
    align <- GRangesList(GRanges("chr1", IRanges(5, 10)),
                         GRanges("chr1", IRanges(c(5, 15), width=6)),
                         GRanges("chr1", IRanges(5, 10)))
    names(align) <- letters[seq_along(align)]

    strand(x) <- "-"
    strand(align) <- "+"
    ans <- pmapToTranscripts(x, align, FALSE)
    checkTrue(length(x) == length(x))
    checkTrue(all(width(ans) == 0L))
    checkIdentical(names(ans), LETTERS[seq_along(align)])

    strand(align[[2]][1]) <- "-"
    checkException(pmapToTranscripts(x, align, FALSE), silent=TRUE) 

    strand(align) <- "+"
    strand(x) <- "+"
    ans <- pmapToTranscripts(x, align, FALSE) 
    checkIdentical(width(ans), c(1L, 1L, 0L))
    checkIdentical(start(ans), c(2L, 8L, 0L))

    strand(align) <- "-"
    strand(x) <- "-"
    ans <- pmapToTranscripts(x, align, FALSE) 
    checkIdentical(width(ans), c(1L, 1L, 0L))
    checkIdentical(start(ans), c(5L, 5L, 0L))
    checkIdentical(names(ans), names(x))
    checkIdentical(seqlevels(ans), c("a", "b", "UNMAPPED"))
    checkIdentical(as.character(strand(ans)), c("-", "-", "*"))

    ## out of bounds
    x <- GRanges("chr1", IRanges(rep(40, 4), width=11))
    align <- GRanges("chr1", IRanges(c(1, 35, 45, 55), width=11))
    ans <- pmapToTranscripts(x, align) 
    checkIdentical(seqnames(ans), Rle(as.factor("UNMAPPED"), 4)) 
}

test_pmapFromTranscripts <- function()
{
    ## empty
    ans1 <- pmapFromTranscripts(GRanges(), GRanges())
    ans2 <- pmapFromTranscripts(GRanges(), GRangesList())
    checkTrue(length(ans1) == 0L)
    checkTrue(length(ans2) == 0L)
    checkIdentical(names(mcols(ans1)), character(0))
    checkIdentical(names(mcols(ans2)), character(0))

    ## length
    x <- GRanges("chr1", IRanges(1, width=1))
    align <- GRanges("chr1", IRanges(1:5, width=1))
    checkException(pmapFromTranscripts(x, align), silent=TRUE)
    align <-  GRangesList(align, align)
    checkException(pmapFromTranscripts(x, align), silent=TRUE)

    ## strand
    x <- GRanges("chr1", IRanges(c(1, 50, 150), width=1, names=LETTERS[1:3]))
    gr <- GRanges("chr1", IRanges(c(100, 300), width=100))
    align <- GRangesList(gr, gr, gr)

    strand(x) <- "-"
    strand(align) <- "+"
    ans <- pmapFromTranscripts(x, align)
    checkTrue(length(x) == length(x))
    checkTrue(all(width(ans) == 0L))

    strand(align[[2]][1]) <- "-"
    checkException(pmapFromTranscripts(x, align), silent=TRUE) 

    strand(align) <- "+"
    strand(x) <- "+"
    ans <- pmapFromTranscripts(x, align) 
    checkIdentical(start(ans), c(100L, 149L, 349L))

    strand(align) <- "-"
    strand(x) <- "-"
    ans <- pmapFromTranscripts(x, align) 
    checkIdentical(start(ans), c(399L, 350L, 150L))
    checkIdentical(names(ans), names(x))
    checkIdentical(seqlevels(ans), "chr1")

    ## out of bounds
    x <- GRanges("chr1", IRanges(rep(40, 3), width=11))
    align <- GRanges("chr1", IRanges(c(35, 45, 55), width=20))
    ans <- pmapFromTranscripts(x, align, TRUE) 
    checkIdentical(seqnames(ans), Rle(as.factor("UNMAPPED"), 3)) 
}
