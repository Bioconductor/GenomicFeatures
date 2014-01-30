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
        stop("'strand' must be a vector, factor, Rle, or list-like object")
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
                 "'transcripts' must be a RangesList object")
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
        stop("when 'x' is a BSgenome object, ",
             "'transcripts' must be a GRangesList or\n  TranscriptDb object")
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

.extractTranscriptSeqsFromOneBSgenomeSeq <-
    function(seqname, x, transcripts)
{
    seqlevels(transcripts, force=TRUE) <- seqname
    strand <- strand(transcripts)
    transcripts <- ranges(transcripts)
    ## Load the sequence.
    if (seqname %in% seqnames(x)) {
        x <- x[[seqname]]
        masks(x) <- NULL
    } else {
        regex <- paste0("^", seqname, "$")
        x <- getSeq(x, regex, as.character=FALSE)
    }
    extractTranscriptSeqs(x, transcripts, strand=strand)
}

setMethod("extractTranscriptSeqs", "BSgenome",
    function(x, transcripts)
    {
        transcripts <- .normarg_transcripts(transcripts)
        seqlevels(transcripts) <- seqlevelsInUse(transcripts)
        ## 'seqnames0' is just an ordinary factor (not Rle) parallel to
        ## 'transcripts'.
        seqnames0 <- unlist(runValue(seqnames(transcripts)), use.names=FALSE)
        dnaset_list <- lapply(levels(seqnames0),
                              .extractTranscriptSeqsFromOneBSgenomeSeq,
                              x, transcripts)
        ans <- unsplit_list_of_XVectorList("DNAStringSet",
                                           dnaset_list,
                                           seqnames0)
        names(ans) <- names(transcripts)
        ans 
    }
)

