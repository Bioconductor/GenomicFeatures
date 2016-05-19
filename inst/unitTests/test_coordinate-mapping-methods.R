library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
cdsbytx <- cdsBy(txdb, "tx", use.names=TRUE)[1:3]
quiet <- suppressWarnings

### ----------------------------------------------------------------------
### mapFromTranscripts(), pmapFromTranscripts()

test_mapFromTranscripts <- function()
{
    ## empty
    ans1 <- mapFromTranscripts(GRanges(), GRanges())
    ans2 <- mapFromTranscripts(GRanges(), GRangesList())
    checkIdentical(ans1, ans2)
    checkIdentical(names(mcols(ans1)), c("xHits", "transcriptsHits"))

    ## 'transcripts' must have names
    x <- GRanges("tx1", IRanges(1, 5), strand="+")
    gr <- GRanges("chr1", IRanges(c(1, 20), width=10), strand="+")
    checkException(mapFromTranscripts(x, gr), silent=TRUE) 
    checkException(mapFromTranscripts(x, GRangesList(gr, gr)), silent=TRUE) 
}

test_mapFromTranscripts_strand <- function()
{
    x <- GRanges("tx1", IRanges(1, 5), strand="+")
    gr <- GRanges("chr1", IRanges(c(1, 20), width=10), strand="+")
    names(gr) <- c("tx1", "tx1")
    align = GRangesList(gr, gr, gr)
    names(align) <- c("tx1", "bar", "tx1") 

    ## "+' strand 
    ans <- mapFromTranscripts(x, align)  ## GRangesList
    checkTrue(all(width(ans) == 5L))
    checkIdentical(mcols(ans)$transcriptsHits, c(1L, 3L))
    checkIdentical(ranges(ans), IRanges(c(1, 1), c(5, 5)))
    checkIdentical(unique(as.character(seqnames(ans))), "chr1")
    ignore <- mapFromTranscripts(x, align, TRUE)
    checkIdentical(ranges(ans), ranges(ignore))
    checkTrue(unique(as.character(strand(ignore))) == "*")
    ans <- mapFromTranscripts(x, gr)  ## GRanges
    checkTrue(all(width(ans) == 5L))
    checkIdentical(mcols(ans)$transcriptsHits, c(1L, 2L))
    checkIdentical(ranges(ans), IRanges(c(1, 20), c(5, 24)))
    checkIdentical(unique(as.character(seqnames(ans))), "chr1")
    ignore <- mapFromTranscripts(x, gr, TRUE)
    checkIdentical(ranges(ans), ranges(ignore))
    checkTrue(unique(as.character(strand(ignore))) == "*")


    ## '-' strand
    xx <- GRanges("tx1", IRanges(12, 14), strand="-")
    gg <- GRanges("chr1", IRanges(c(1, 20), width=10), strand="-")
    aa = GRangesList(gg)
    names(aa) <- "tx1" 
    ans <- mapFromTranscripts(xx, aa)
    checkIdentical(ranges(ans), IRanges(7, 9))

    ## invalid mixed strand 
    strand(align[[1]][1]) <- "-"
    checkException(mapFromTranscripts(x, align, FALSE), silent=TRUE) 

    ## valid mixed strand:
    strand(align[[1]]) <- "-"
    ans <- mapFromTranscripts(x, align) ## GRangesList
    checkTrue(all(width(ans) == 5L))
    checkIdentical(mcols(ans)$transcriptsHits, 3L)
    checkIdentical(ranges(ans), IRanges(1, 5))
    checkIdentical(unique(as.character(seqnames(ans))), "chr1")
    ignore <- mapFromTranscripts(x, align, TRUE)
    checkIdentical(mcols(ignore)$transcriptsHits, c(1L, 3L))
    checkTrue(unique(as.character(strand(ignore))) == "*")
    strand(x) <- "-"
    ans <- mapFromTranscripts(x, align)
    checkIdentical(ranges(ans), IRanges(25, end=29))
    strand(gr[2]) <- "-"
    strand(x) <- "+"
    ans <- mapFromTranscripts(x, gr)  ## GRanges
    checkTrue(all(width(ans) == 5L))
    checkIdentical(mcols(ans)$transcriptsHits, 1L)
    checkIdentical(ranges(ans), IRanges(1, 5))
    checkIdentical(unique(as.character(seqnames(ans))), "chr1")
    ignore <- mapFromTranscripts(x, gr, TRUE)
    checkIdentical(ranges(ignore), IRanges(c(1, 20), c(5, 24)))
    checkTrue(unique(as.character(strand(ignore))) == "*")
    strand(x) <- "-"
    ans <- mapFromTranscripts(x, align)
    checkIdentical(ranges(ans), IRanges(25, end=29))
}

test_mapFromTranscripts_order_in_GRangesList <- function()
{
    x <- GRanges("tx1", IRanges(1, 5), strand="+")
    gr1 <- GRanges("chr1", IRanges(c(1, 20), end=c(10, 30)), strand="+")
    gr2 <- GRanges("chr1", IRanges(c(20, 1), end=c(30, 10)), strand="+")
    align <- GRangesList(gr1, gr2); names(align) <- c("tx1", "tx1")

    ans <- mapFromTranscripts(x, align) 
    checkIdentical(ranges(ans), IRanges(c(1, 1), c(5, 5)))
    checkIdentical(mcols(ans)$xHits, c(1L, 1L))
    checkIdentical(mcols(ans)$transcriptsHits, c(1L, 2L))

    strand(x) <- strand(align) <- "-"
    ans <- mapFromTranscripts(x, align) 
    checkIdentical(ranges(ans), IRanges(c(26, 26), c(30, 30)))
    checkIdentical(mcols(ans)$xHits, c(1L, 1L))
    checkIdentical(mcols(ans)$transcriptsHits, c(1L, 2L))
}

test_pmapFromTranscripts <- function()
{
    ## empty
    ans1 <- pmapFromTranscripts(IRanges(), GRanges())
    ans2 <- pmapFromTranscripts(IRanges(), GRangesList())
    ans3 <- pmapFromTranscripts(GRanges(), GRanges())
    ans4 <- pmapFromTranscripts(GRanges(), GRangesList())
    checkTrue(length(ans1) == length(ans2))
    checkTrue(length(ans3) == length(ans4))
    checkTrue(length(ans1) == 0L)
    checkTrue(is(ans1, "GRanges"))
    checkTrue(is(ans2, "GRangesList"))
    checkTrue(is(ans3, "GRanges"))
    checkTrue(is(ans4, "GRangesList"))
}

test_pmapFromTranscripts_recycling <- function()
{
    x <- GRanges("tx1", IRanges(c(1, 1, 1), width=5), strand="+")
    gr <- GRanges("chr1", IRanges(c(10, 10, 10), width=10), strand="+")

    checkException(pmapFromTranscripts(x[1:2], gr), silent=TRUE)
    checkException(pmapFromTranscripts(x, gr[1:2]), silent=TRUE)

    ## 'transcripts' as GRanges
    ans1 <- pmapFromTranscripts(x[1], gr)
    ans2 <- pmapFromTranscripts(x, gr[1])
    checkIdentical(ans1, ans2)
    ans1 <- pmapFromTranscripts(ranges(x[1]), gr)
    ans2 <- pmapFromTranscripts(ranges(x), gr[1])
    checkIdentical(ans1, ans2)
    checkIdentical(names(mcols(ans1)), character(0))

    ## 'transcripts' as GRangesList
    grl <- split(gr, seq_along(gr))
    ans1 <- pmapFromTranscripts(x[1], grl)
    ans2 <- pmapFromTranscripts(x, grl[1])
    checkIdentical(unname(ans1), unname(ans2))
    ans1 <- pmapFromTranscripts(ranges(x[1]), grl)
    ans2 <- pmapFromTranscripts(ranges(x), grl[1])
    checkIdentical(unname(ans1), unname(ans2))
}

test_pmapFromTranscripts_strand <- function()
{
    x <- GRanges("chr1", IRanges(1, 5))
    align <- GRangesList(
        GRanges("chr1", IRanges(c(201, 101), c(220, 120)), strand="-"),
        GRanges("chr1", IRanges(501, 535), strand="+"))

    ans <- pmapFromTranscripts(ranges(x), align)
    checkIdentical(width(ans), IntegerList(c(5, 0), 5))
    checkIdentical(start(ans), IntegerList(c(216, 101), 501))

    ## 'x' is "*"
    ans <- pmapFromTranscripts(x, align)
    checkIdentical(width(ans), IntegerList(c(5, 0), 5))
    checkIdentical(start(ans), IntegerList(c(216, 101), 501))

    ## 'x' is "+"
    strand(x) <- "+"
    ans <- pmapFromTranscripts(x, align)
    checkIdentical(width(ans), IntegerList(c(0, 0), 5))
    checkIdentical(start(ans), IntegerList(c(201, 101), 501))

    ## 'x' is "-"
    strand(x) <- "-"
    ans <- pmapFromTranscripts(x, align)
    checkIdentical(width(ans), IntegerList(c(5, 0), 0))
    checkIdentical(start(ans), IntegerList(c(216, 101), 501))

    ## ignore.strand
    ans <- pmapFromTranscripts(x, align, TRUE)
    checkIdentical(width(ans), IntegerList(c(0, 5), 5))
    checkIdentical(start(ans), IntegerList(c(201, 101), 501))
    checkIdentical(as.character(runValue(strand(ans))), c("*", "*")) 

    ## order of ranges in return GRangesList 
    strand(align) <- strand(x) <- "+"
    ans <- pmapFromTranscripts(x, align)

    ## invalid mixed strand
    x <- GRanges("chr1", IRanges(1, 5), strand="+")
    align <- GRangesList(
        GRanges("chr1", IRanges(c(1, 1), width=5), strand=c("+", "-")))
    checkException(pmapFromTranscripts(x, align), silent=TRUE)
    checkException(pmapFromTranscripts(ranges(x), align), silent=TRUE)
}

test_pmapFromTranscripts_UNMAPPED <- function()
{
    ## 'transcripts' is GRanges
    x <- GRanges("chr1", IRanges(40, 50))
    align <- GRanges("chr1", IRanges(35, width=5))
    ans <- pmapFromTranscripts(x, align, TRUE) 
    checkIdentical(as.character(seqnames(ans)), "UNMAPPED")
    checkTrue(width(ans) == 0L)
}

### ----------------------------------------------------------------------
### mapToTranscripts(), pmapToTranscripts()

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

    ans <- mapToTranscripts(x, cdsbytx, ignore.strand=FALSE)
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

    x <- GRanges("chr3", IRanges(9, 9), strand="+")
    transcripts <- GRangesList(tx1=GRanges("chr3", IRanges(3, 10), strand="-"))
    ans <- mapToTranscripts(x, transcripts)
    checkTrue(length(ans) == 0L)
    ans <- mapToTranscripts(x, transcripts, TRUE)
    checkTrue(start(ans) == 7L)
    checkTrue(as.character(strand(ans)) == "*")

    x <- GRanges("1", IRanges(248, width=1))
    transcripts <- GRangesList(
        foo = (GRanges("1", IRanges(c(101,201), width=50), strand="-")))
    ans <- mapToTranscripts(x, transcripts, TRUE)
    checkIdentical(ranges(ans), IRanges(98, 98))
    checkTrue(as.character(strand(ans)) == "*")
    ans <- mapToTranscripts(x, transcripts, FALSE)
    checkIdentical(ranges(ans), IRanges(3, 3))
    checkTrue(as.character(strand(ans)) == "-")

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

    ## methods  GR,GR  GR,GRL  
    ## recycling 
    x <- GRanges("chr1", IRanges(1, width=1))
    y <- GRanges("chr1", IRanges(1:5, width=1))
    ans <- pmapToTranscripts(x, GRangesList("tx1"=y, "tx2"=y))
    checkIdentical(as.character(seqnames(ans)), c("tx1", "tx2"))
    ans <- pmapToTranscripts(c(x, x), GRangesList("tx1"=y))
    checkIdentical(as.character(seqnames(ans)), c("tx1", "tx1"))

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

test_mapToTranscripts_intronJunctions <- function()
{
    ## 2 out of bounds, intron, true hit
    x <- GRanges("chr1", IRanges(start=c(1, 6, 8, 15), width=1))
    align <- GRangesList(
        "A"=GRanges("chr1", IRanges(start=c(5, 10), end=c(6, 12))))
    ans <- mapToTranscripts(x, align, intronJunctions=TRUE)
    checkTrue(length(ans) == 2)
    checkEquals(mcols(ans)$xHits, c(2, 3))
    checkEquals(width(ans), c(1, 0))
    checkEquals(start(ans), c(2, 3))
    checkEquals(mcols(ans)$intronic, c(FALSE, TRUE))

    ans <- mapToTranscripts(rev(x), align, intronJunctions=TRUE)
    checkEquals(width(ans), c(0, 1))
    checkEquals(mcols(ans)$intronic, c(TRUE, FALSE))

    ## span intron and exon
    x <- GRanges("chr1", IRanges(start=6, end=7))
    ans <- mapToTranscripts(x, align, intronJunctions=TRUE)
    checkTrue(length(ans) == 0)

    ## no non-hits
    x <- GRanges("chr1", IRanges(start=5, end=6))
    ans <- mapToTranscripts(x, align, intronJunctions=TRUE)
    checkTrue(ranges(ans) == IRanges(1, 2))
    checkEquals(mcols(ans)$intronic, logical(1))

    ## negative strand
    x <- GRanges("chr1", IRanges(start=c(1, 6, 8, 15), width=1))
    align <- GRangesList(
        "A"=GRanges("chr1", IRanges(start=c(5, 10), end=c(6, 12))))
    strand(x) <- strand(align) <- "-"
    ans <- mapToTranscripts(x, align, intronJunctions=TRUE)
    checkTrue(length(ans) == 2)
    checkEquals(width(ans), c(1, 0))
    checkEquals(start(ans), c(4, 4))
}
