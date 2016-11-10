### =========================================================================
### extractTranscriptSeqs() and related tools
### -------------------------------------------------------------------------


.unlist_strand <- function(strand, transcripts)
{
    if (is.list(strand) || is(strand, "List")) {
        ## 'strand' is a list-like object.
        if (!identical(unname(elementNROWS(strand)),
                       unname(elementNROWS(transcripts))))
            stop(wmsg("when 'strand' is a list-like object, it must have ",
                      "the same \"shape\" as 'transcripts' (i.e. same length ",
                      "plus 'strand[[i]]' must have the same length as ",
                      "'transcripts[[i]]' for all 'i')"))
        return(strand(unlist(strand, use.names=FALSE)))
    }
    if (!(is.vector(strand) || is.factor(strand) || is(strand, "Rle")))
        stop(wmsg("'strand' must be an atomic vector, a factor, ",
                  "an Rle object, or a list-like object"))
    strand <- strand(strand)
    strand <- S4Vectors:::V_recycle(strand, transcripts,
                                    "strand", "transcripts")
    rep.int(strand, elementNROWS(transcripts))
}

setGeneric("extractTranscriptSeqs", signature="x",
    function(x, transcripts, ...) standardGeneric("extractTranscriptSeqs")
)

setMethod("extractTranscriptSeqs", "DNAString",
    function(x, transcripts, strand="+")
    {
        if (!is(transcripts, "RangesList"))
            stop(wmsg("when 'x' is a DNAString object, ",
                      "'transcripts' must be an RangesList object"))
        unlisted_strand <- .unlist_strand(strand, transcripts)
        if (!all(unlisted_strand %in% c("+", "-")))
            stop(wmsg("'strand' can only contain \"+\" and/or \"-\" values. ",
                      "\"*\" is not allowed."))
        idx <- which(unlisted_strand == "-")
        exons <- extractList(x, unlist(transcripts, use.names=FALSE))
        exons[idx] <- reverseComplement(exons[idx])
        unstrsplit(relist(exons, transcripts))
    }
)

### Check for transcripts that have exons located on more than one
### chromosome.
.check_exon_chrom <- function(tx1)
{
    run_lens <- runLength(seqnames(tx1))
    idx <- which(elementNROWS(run_lens) != 1L)
    if (length(idx) == 0L)
        return()
    tx1_names <- names(tx1)
    if (is.null(tx1_names)) {
        some_in1string <- ""
    } else {
        some_idx <- head(idx, n=2L)
        some_names <- tx1_names[some_idx]
        some_in1string <- paste0(some_names, collapse=", ")
        if (length(idx) > length(some_idx))
            some_in1string <- paste0("e.g. ", some_in1string, ", etc...")
        some_in1string <- paste0(" (", some_in1string, ")")
    }
    stop(wmsg("Some transcripts", some_in1string, " have exons located on ",
              "more than one chromosome. This is not supported yet."))
}

### Check the "exon_rank" inner metadata column if present. When 'transcripts'
### contains CDSs (instead of exons) grouped by transcript, some of the lowest
### or/and highest exon ranks can be missing.
.check_exon_rank <- function(tx1)
{
    exon_rank <- mcols(tx1@unlistData)$exon_rank
    if (is.null(exon_rank))
        return()
    if (!is.numeric(exon_rank))
        stop(wmsg("\"exon_rank\" inner metadata column in GRangesList ",
                  "object 'transcripts' is not numeric"))
    if (!is.integer(exon_rank)) {
        warning(wmsg("\"exon_rank\" inner metadata column in GRangesList ",
                     "object 'transcripts' is not integer"))
        exon_rank <- as.integer(exon_rank)
    }
    if (any(is.na(exon_rank)))
        stop(wmsg("\"exon_rank\" inner metadata column in GRangesList ",
                  "object 'transcripts' contains NAs"))

    partitioning <- PartitioningByEnd(tx1)
    ## The 2 lines below are equivalent to:
    ##   tmp <- relist(exon_rank, partitioning)
    ##   min_rank <- min(tmp)
    ## but much faster!
    v <- Views(exon_rank, partitioning)
    min_rank <- viewMins(v)
    if (any(min_rank < 1L))
        stop(wmsg("\"exon_rank\" inner metadata column in GRangesList ",
                  "object 'transcripts' contains ranks < 1"))
    tx1_eltNROWS <- elementNROWS(partitioning)
    target <- S4Vectors:::fancy_mseq(tx1_eltNROWS,
                                     offset=min_rank - 1L)
    if (!identical(target, unname(exon_rank)))
        stop(wmsg("\"exon_rank\" inner metadata column in GRangesList ",
                  "object 'transcripts' does not contain increasing ",
                  "consecutive ranks for some transcripts"))
}

### TODO: Incorporate this fast path to "unlist" method for XStringSet objects.
.fast_XStringSet_unlist <- function(x)
{
    # Disabling the fast path for now. Until I understand why using it
    # causes extractTranscriptSeqs(Hsapiens, TxDb.Hsapiens.UCSC.hg18.knownGene)
    # to use more memory (319.7 Mb) than when NOT using it (288.9 Mb).
if (FALSE) {
    x_len <- length(x)
    if (x_len != 0L && length(x@pool) == 1L) {
        x_ranges <- x@ranges
        x_start <- start(x_ranges)
        x_end <- end(x_ranges)
        if (identical(x_end[-x_len] + 1L, x_start[-1L])) {
            ## The ranges are adjacent. We can unlist() without copying
            ## the sequence data!
            cat("using fast path (", x_len, ") ...\n")
            ans_class <- elementType(x)
            ans_shared <- x@pool[[1L]]
            ans_offset <- x_start[1L] - 1L
            ans_length <- x_end[x_len] - ans_offset
            ans <- new2(ans_class, shared=ans_shared,
                                   offset=ans_offset, 
                                   length=ans_length,
                                   check = FALSE)
            return(ans)
        }
    }
}
    unlist(x, use.names=FALSE)
}

.extract_and_combine <- function(x, seqname, ranges)
{
    seqs <- getSeq(x, GRanges(seqname, ranges))
    ## For "getSeq" methods (like the method for GmapGenome objects) that
    ## return a character vector.
    if (is.character(seqs))
        seqs <- DNAStringSet(seqs)
    .fast_XStringSet_unlist(seqs)
}

.extractTranscriptSeqsFromOneSeq <- function(seqlevel, x, transcripts)
{
    seqlevels(transcripts, pruning.mode="coarse") <- seqlevel
    strand <- strand(transcripts)
    transcripts <- ranges(transcripts)
    if (seqlevel %in% seqlevels(x)) {
        ## We try to load the less stuff possible i.e. only the nucleotides
        ## that participate in at least one exon.
        exons <- unlist(transcripts, use.names=FALSE)
        ranges_to_load <- reduce(exons, with.inframe.attrib=TRUE)
        x <- .extract_and_combine(x, seqlevel, ranges_to_load)
        exons <- attr(ranges_to_load, "inframe")
        transcripts <- relist(exons, transcripts)
    } else {
        ## Why do we need this?
        regex <- paste0("^", seqlevel, "$")
        x <- getSeq(x, regex)
    }
    extractTranscriptSeqs(x, transcripts, strand=strand)
}

setMethod("extractTranscriptSeqs", "ANY",
    function(x, transcripts, ...)
    {
        if (is(transcripts, "GRangesList")) {
            if (length(list(...)) != 0L)
                stop(wmsg("additional arguments are allowed only when ",
                          "'transcripts' is not a GRangesList object"))
        } else {
            transcripts <- try(exonsBy(transcripts, by="tx", ...),
                               silent=TRUE)
            if (is(transcripts, "try-error"))
                stop(wmsg("failed to extract the exon ranges ",
                          "from 'transcripts' ",
                          "with exonsBy(transcripts, by=\"tx\", ...)"))
        }
        idx1 <- which(elementNROWS(transcripts) != 0L)
        tx1 <- transcripts[idx1]
        .check_exon_chrom(tx1)
        .check_exon_rank(tx1)

        seqlevels(tx1) <- seqlevelsInUse(tx1)
        ## 'seqnames1' is just an ordinary factor (not Rle) parallel to 'tx1'.
        seqnames1 <- unlist(runValue(seqnames(tx1)), use.names=FALSE)
        dnaset_list <- lapply(levels(seqnames1),
                              .extractTranscriptSeqsFromOneSeq, x, tx1)
        ans <- rep.int(DNAStringSet(""), length(transcripts))
        names(ans) <- names(transcripts)
        ans[idx1] <- unsplit_list_of_XVectorList("DNAStringSet",
                                                 dnaset_list,
                                                 seqnames1)
        ans 
    }
)

