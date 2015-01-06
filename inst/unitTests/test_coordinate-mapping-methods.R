library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
cdsbytx <- cdsBy(txdb, "tx")[1:3]

test_mapToTranscript <- function()
{
    ## reverse = FALSE

    ## empty
    ans1 <- mapToTranscript(GRanges(), GRanges())
    ans2 <- mapToTranscript(GRanges(), GRangesList())
    checkTrue(length(ans1) == 0L)
    checkTrue(length(ans2) == 0L)
    checkIdentical(names(mcols(ans1)), c("xHits", "transcriptHits"))
    checkIdentical(names(mcols(ans2)), c("xHits", "transcriptHits"))

    ## strand
    x <- GRanges("chr2L", IRanges(c(7500, 8400, 9000), 
                 width=200, names=LETTERS[1:3]), strand="+")
    align <- cdsbytx

    ans <- mapToTranscript(x, cdsbytx, FALSE, FALSE)
    checkTrue(length(ans) == 3L)
    checkIdentical(names(ans), c("B", "B", "C"))
    checkTrue(all(width(ans) == 200L))
    checkIdentical(mcols(ans)$xHits, c(2L, 2L, 3L))
    checkIdentical(mcols(ans)$transcriptHits, c(1L, 3L, 2L))

    strand(align[[2]][1]) <- "-"
    checkException(mapToTranscript(x, align, ignore.strand=FALSE), silent=TRUE) 

    strand(align[[2]]) <- "-"
    ans <- mapToTranscript(x, align, ignore.strand=FALSE) 
    checkIdentical(mcols(ans)$xHits, c(2L, 2L))
    checkIdentical(mcols(ans)$transcriptHits, c(1L, 3L))

    ## TODO: seqnames

    ## reverse = TRUE

    ## strand
    x <- GRanges("chr1", IRanges(rep(10, 5), width=2), strand="+")
    names(x) <- c("A", "B", "C", "foo", "A")
    gr <- GRanges("chr1", IRanges(c(10, 100), width=50), strand="+")
    align = GRangesList(gr, gr, gr, gr)
    names(align) <- c("C", "bar", "C", "A") 

    ans <- mapToTranscript(x, align, TRUE, TRUE)
    checkTrue(length(ans) == 4L)
    checkIdentical(names(ans), c("A", "C", "C", "A"))
    checkTrue(all(width(ans) == 2L))
    checkIdentical(mcols(ans)$xHits, c(1L, 3L, 3L, 5L))
    checkIdentical(mcols(ans)$transcriptHits, c(4L, 1L, 3L, 4L))
    ans <- mapToTranscript(x, unlist(align), TRUE, TRUE)
    checkTrue(length(ans) == 8L)
    checkIdentical(names(ans), c("A", "A", rep("C", 4), "A", "A"))
    checkIdentical(mcols(ans)$transcriptHits, c(7L, 8L, 1L, 2L, 5L, 6L, 7L, 8L))

    strand(align[[1]][1]) <- "-"
    checkException(mapToTranscript(x, align, TRUE, FALSE), silent=TRUE) 

    strand(align[[1]]) <- "-"
    ans <- mapToTranscript(x, align, TRUE, FALSE) 
    checkIdentical(mcols(ans)$transcriptHits, c(4L, 3L, 4L))
    ans <- mapToTranscript(x, align, TRUE, TRUE) 
    checkIdentical(mcols(ans)$transcriptHits, c(4L, 1L, 3L, 4L))

    ## TODO: seqnames
}

test_mapToTranscript_range_order <- function()
{
    x <- GRanges("chrA", IRanges(43522349, width=1), strand="+")
    align1 <- GRangesList(GRanges("chrA", 
        IRanges(c(43522244, 43528406),
                c(43524145, 43528644)), strand="+"))
    align2 <- GRangesList(GRanges("chrA", 
        IRanges(c(43528406, 43522244),
                c(43528644, 43524145)), strand="+"))

    ## "+" strand
    ## smallest range first 
    ans <- mapToTranscript(x, align1, FALSE, FALSE)
    checkTrue(start(ans) == 106L)
    ## largest range first 
    ans <- mapToTranscript(x, align2, FALSE, FALSE)
    checkTrue(start(ans) == 106L)

    ## "-" strand
    strand(x) <- "-"
    strand(align1) <- "-"
    strand(align2) <- "-"
    ## smallest range first
    ans <- mapToTranscript(x, align1, FALSE, FALSE)
    checkTrue(start(ans) == 2036L)
    ## largest range first
    ans <- mapToTranscript(x, align2, FALSE, FALSE)
    checkTrue(start(ans) == 2036L)
}

test_pmapToTranscript <- function()
{
    ## reverse = FALSE

    ## empty
    ans1 <- pmapToTranscript(GRanges(), GRanges())
    ans2 <- pmapToTranscript(GRanges(), GRangesList())
    checkTrue(length(ans1) == 0L)
    checkTrue(length(ans2) == 0L)
    checkIdentical(names(mcols(ans1)), character(0))
    checkIdentical(names(mcols(ans2)), character(0))

    ## length
    x <- GRanges("chr1", IRanges(1, width=1))
    y <- GRanges("chr1", IRanges(1:5, width=1))
    checkException(pmapToTranscript(x, y), silent=TRUE)
    checkException(pmapToTranscript(x, GRangesList(y, y)), silent=TRUE)

    ## strand
    x <- GRanges("chr1", IRanges(c(6, 16, 1), width=1, names=LETTERS[1:3]))
    align <- GRangesList(GRanges("chr1", IRanges(5, 10)),
                         GRanges("chr1", IRanges(c(5, 15), width=6)),
                         GRanges("chr1", IRanges(5, 10)))

    strand(x) <- "-"
    strand(align) <- "+"
    ans <- pmapToTranscript(x, align, FALSE, FALSE)
    checkTrue(length(x) == length(x))
    checkTrue(all(width(ans) == 0L))
    checkIdentical(names(ans), LETTERS[1:3])

    strand(align[[2]][1]) <- "-"
    checkException(pmapToTranscript(x, align, FALSE, FALSE), silent=TRUE) 

    strand(align) <- "+"
    strand(x) <- "+"
    ans <- pmapToTranscript(x, align, FALSE, FALSE) 
    checkIdentical(width(ans), c(1L, 1L, 0L))
    checkIdentical(start(ans), c(2L, 8L, 1L))

    strand(align) <- "-"
    strand(x) <- "-"
    ans <- pmapToTranscript(x, align, FALSE, FALSE) 
    checkIdentical(width(ans), c(1L, 1L, 0L))
    checkIdentical(start(ans), c(5L, 5L, 1L))

    ## out of bounds
    x <- GRanges("chr1", IRanges(rep(40, 4), width=11))
    align <- GRanges("chr1", IRanges(c(1, 35, 45, 55), width=11))
    ans <- pmapToTranscript(x, align) 
    checkIdentical(seqnames(ans), Rle(as.factor("unmapped"), 4)) 

    ## reverse = TRUE

    ## empty
    ans1 <- pmapToTranscript(GRanges(), GRanges())
    ans2 <- pmapToTranscript(GRanges(), GRangesList())
    checkTrue(length(ans1) == 0L)
    checkTrue(length(ans2) == 0L)
    checkIdentical(names(mcols(ans1)), character(0))
    checkIdentical(names(mcols(ans2)), character(0))

    ## length
    x <- GRanges("chr1", IRanges(1, width=1))
    align <- GRanges("chr1", IRanges(1:5, width=1))
    checkException(pmapToTranscript(x, align), silent=TRUE)
    checkException(pmapToTranscript(x, GRangesList(align, align)), silent=TRUE)

    ## strand
    x <- GRanges("chr1", IRanges(c(1, 50, 150), width=1, names=LETTERS[1:3]))
    gr <- GRanges("chr1", IRanges(c(100, 300), width=100))
    align <- GRangesList(gr, gr, gr)

    strand(x) <- "-"
    strand(align) <- "+"
    ans <- pmapToTranscript(x, align, TRUE, FALSE)
    checkTrue(length(x) == length(x))
    checkTrue(all(width(ans) == 0L))
    checkIdentical(names(ans), LETTERS[1:3])

    strand(align[[2]][1]) <- "-"
    checkException(pmapToTranscript(x, align, TRUE, FALSE), silent=TRUE) 

    strand(align) <- "+"
    strand(x) <- "+"
    ans <- pmapToTranscript(x, align, TRUE, FALSE) 
    checkIdentical(start(ans), c(100L, 149L, 349L))

    strand(align) <- "-"
    strand(x) <- "-"
    ans <- pmapToTranscript(x, align, TRUE, FALSE) 
    checkIdentical(start(ans), c(399L, 350L, 150L))

    ## out of bounds
    x <- GRanges("chr1", IRanges(rep(40, 3), width=11))
    align <- GRanges("chr1", IRanges(c(35, 45, 55), width=20))
    ans <- pmapToTranscript(x, align, TRUE, TRUE) 
    checkIdentical(seqnames(ans), Rle(as.factor("unmapped"), 3)) 
}
