### =========================================================================
### extractTranscriptSeqs() and related tools
### -------------------------------------------------------------------------


.unlist_strand <- function(strand, transcripts)
{
    if (is.list(strand) || is(strand, "List")) {
        ## 'strand' is a list-like object.
        if (!identical(unname(elementLengths(strand)),
                       unname(elementLengths(transcripts))))
            stop("when 'strand' is a list-like object, it must have ",
                 "the same \"shape\" as\n  'transcripts' (i.e. same length ",
                 "plus 'strand[[i]]' must have the same\n  length as ",
                 "'transcripts[[i]]' for all 'i')")
        return(strand(unlist(strand, use.names=FALSE)))
    }
    if (!(is.vector(strand) || is.factor(strand) || is(strand, "Rle")))
        stop("'strand' must be an atomic vector, a factor, an Rle object, ",
             "or a list-like object")
    strand <- strand(strand)
    strand <- Biostrings:::.V_recycle(strand, transcripts,
                                      "strand", "length of 'transcripts'")
    rep.int(strand, elementLengths(transcripts))
}

setGeneric("extractTranscriptSeqs", signature="x",
    function(x, transcripts, ...) standardGeneric("extractTranscriptSeqs")
)

setMethod("extractTranscriptSeqs", "DNAString",
    function(x, transcripts, strand="+")
    {
        if (!is(transcripts, "RangesList"))
            stop("when 'x' is a DNAString object, ",
                 "'transcripts' must be an RangesList object")
        unlisted_strand <- .unlist_strand(strand, transcripts)
        if (!all(unlisted_strand %in% c("+", "-")))
            stop("'strand' can only contain \"+\" and/or \"-\" values. ",
                 "\"*\" is not allowed.")
        idx <- which(unlisted_strand == "-")
        exons <- extractList(x, unlist(transcripts, use.names=FALSE))
        exons[idx] <- reverseComplement(exons[idx])
        unstrsplit(relist(exons, transcripts))
    }
)

## When 'transcripts' contains CDSs (instead of exons) grouped by transcript,
## some of the lowest or/and highest exon ranks can be missing.
.check_exon_rank <- function(exon_rank, partitioning)
{
    if (!is.numeric(exon_rank))
        stop("\"exon_rank\" column in GRangesList object 'transcripts' ",
             "is not numeric")
    if (!is.integer(exon_rank)) {
        warning("\"exon_rank\" column in GRangesList object 'transcripts' ",
                "is not integer")
        exon_rank <- as.integer(exon_rank)
    }
    ## The 2 lines below are equivalent to:
    ##   tmp <- relist(exon_rank, partitioning)
    ##   min_rank <- min(tmp)
    ## but much faster!
    v <- Views(exon_rank, partitioning)
    min_rank <- viewMins(v)
    if (any(min_rank < 1L))
        stop("\"exon_rank\" column in GRangesList object 'transcripts' ",
             "contains ranks < 1")
    transcripts_eltlens <- elementLengths(partitioning)
    target <- IRanges:::fancy_mseq(transcripts_eltlens,
                                   offset=min_rank - 1L)
    if (!identical(target, unname(exon_rank)))
        stop("\"exon_rank\" column in GRangesList object 'transcripts' ",
             "does not contain\n  increasing consecutive ranks ",
             "for some transcripts")
}

.normarg_transcripts <- function(transcripts)
{
    if (is(transcripts, "TranscriptDb")) {
        transcripts <- exonsBy(transcripts, by="tx", use.names=TRUE)
    } else if (!is(transcripts, "GRangesList")) {
        stop("'transcripts' must be a GRangesList or\n  TranscriptDb object")
    }
    ## Check for transcripts that have exons located on more than one
    ## chromosome.
    run_lens <- runLength(seqnames(transcripts))
    idx <- which(elementLengths(run_lens) != 1L)
    if (length(idx) != 0L) {
        transcripts_names <- names(transcripts)
        if (is.null(transcripts_names)) {
            some_in1string <- ""
        } else {
            some_idx <- head(idx, n=2L)
            some_names <- transcripts_names[some_idx]
            some_in1string <- paste0(some_names, collapse=", ")
            if (length(idx) > length(some_idx))
                some_in1string <- paste0("e.g. ", some_in1string, ", etc...")
            some_in1string <- paste0(" (", some_in1string, ")")
        }
        stop("Some transcripts", some_in1string, " have exons located on ",
             "more than one chromosome.\n  This is not supported yet.")
    }
    ## Check the "rank" metadata column if present.
    exon_rank <- mcols(transcripts@unlistData)$exon_rank
    if (!is.null(exon_rank))
        .check_exon_rank(exon_rank, PartitioningByEnd(transcripts))
    transcripts
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
    .fast_XStringSet_unlist(getSeq(x, GRanges(seqname, ranges)))

.extractTranscriptSeqsFromOneSeq <-
    function(seqname, x, transcripts)
{
    seqlevels(transcripts, force=TRUE) <- seqname
    strand <- strand(transcripts)
    transcripts <- ranges(transcripts)
    if (seqname %in% seqnames(x)) {
        ## We try to load the less stuff possible i.e. only the nucleotides
        ## that participate in at least one exon.
        exons <- unlist(transcripts, use.names=FALSE)
        ranges_to_load <- reduce(exons, with.inframe.attrib=TRUE)
        x <- .extract_and_combine(x, seqname, ranges_to_load)
        exons <- attr(ranges_to_load, "inframe")
        transcripts <- relist(exons, transcripts)
    } else {
        ## Why do we need this?
        regex <- paste0("^", seqname, "$")
        x <- getSeq(x, regex, as.character=FALSE)
    }
    extractTranscriptSeqs(x, transcripts, strand=strand)
}

setMethod("extractTranscriptSeqs", "ANY",
    function(x, transcripts)
    {
        transcripts <- .normarg_transcripts(transcripts)
        seqlevels(transcripts) <- seqlevelsInUse(transcripts)
        ## 'seqnames0' is just an ordinary factor (not Rle) parallel to
        ## 'transcripts'.
        seqnames0 <- unlist(runValue(seqnames(transcripts)), use.names=FALSE)
        dnaset_list <- lapply(levels(seqnames0),
                              ..extractTranscriptSeqsFromOneSeq,
                              x, transcripts)
        ans <- unsplit_list_of_XVectorList("DNAStringSet",
                                           dnaset_list,
                                           seqnames0)
        names(ans) <- names(transcripts)
        ans 
    }
)

