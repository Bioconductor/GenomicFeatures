### =========================================================================
### extractTranscriptCoverage()
### -------------------------------------------------------------------------


### A very efficient way to concatenate groups of successive list elements
### in 'x'.
### 'x' must be a list-like object (typically a CompressedList object).
### 'skeleton' will be immediately replaced with 'PartitioningByEnd(skeleton)'
### so it should be an object that is accepted by the PartitioningByEnd()
### constructor (note that this constructor is a no-op if 'skeleton' is
### already a PartitioningByEnd object). After this transformation, 'skeleton'
### should define a partitioning that is compatible with 'seq_along(x)' (i.e.
### 'relist(seq_along(x), skeleton)' should work).
### Return a list-like object of the same elementType() as 'x' but with the
### length and names of 'skeleton'. Behaves as an endomorphism on a
### CompressedList or PartitioningByEnd object.
### Note that
###     regroup(x, length(x))[[1L]]
### is equivalent to
###     unlist(x, use.names=FALSE)
### TODO:
### - Move this helper to IRanges (R/extractList.R might be a good place) and
###   maybe export and document it?
### - Remove older versions of this helper from GenomicFeatures/R/transcripts.R
###   and SplicingGraphs/R/utils.R, and use IRanges::regroup() instead.
regroup <- function(x, skeleton)
{
    skeleton <- PartitioningByEnd(skeleton)
    x_breakpoints <- end(PartitioningByEnd(x))
    ans_breakpoints <- x_breakpoints[end(skeleton)]
    nleading0s <- length(skeleton) - length(ans_breakpoints)
    if (nleading0s != 0L)
        ans_breakpoints <- c(rep.int(0L, nleading0s), ans_breakpoints)
    ans_partitioning <- PartitioningByEnd(ans_breakpoints,
                                          names=names(skeleton))
    if (is(x, "PartitioningByEnd"))
        return(ans_partitioning)
    relist(unlist(x, use.names=FALSE), ans_partitioning)
}

### The function is named after extractTranscriptSeqs() because extracting
### transcript coverage from a set of aligned reads is an operation that feels
### a lot like extracting transcript sequences from a genome.
### We define it as an ordinary function though (and not as a generic like
### extractTranscriptSeqs()), at least for now.
extractTranscriptCoverage <- function(reads, transcripts)
{
    if (!is(transcripts, "GRangesList")) {
        transcripts <- try(exonsBy(transcripts, by="tx", use.names=TRUE),
                           silent=TRUE)
        if (is(transcripts, "try-error"))
            stop(wmsg("failed to extract the exon ranges ",
                      "from 'transcripts' with ",
                      "exonsBy(transcripts, by=\"tx\", use.names=TRUE)"))
    }
    seqinfo(reads) <- merge(seqinfo(reads), seqinfo(transcripts))

  ## 1) Compute unique exons ('uex').

    ex <- unlist(transcripts, use.names=FALSE)
    ## We could simply do 'uex <- unique(ex)' here but we're going to need
    ## 'sm' and 'is_unique' later to compute the "reverse index" so we compute
    ## them now and use them to extract the unique exons. That way we hash
    ## 'ex' only once (the expensive operation).
    sm <- selfmatch(ex)  # uses a hash table internally
    is_unique <- sm == seq_along(sm)
    uex2ex <- which(is_unique)  # index of unique exons
    uex <- ex[uex2ex]  # unique exons

  ## 2) Compute coverage for each unique exon ('uex_cvg').

    #There doesn't seem to be much benefit in doing this.
    #reads <- subsetByOverlaps(reads, transcripts, ignore.strand=TRUE)
    cvg <- coverage(reads)
    uex_cvg <- cvg[uex]  # parallel to 'uex'

  ## 3) Flip coverage for exons on minus strand.

    ## It feels like this is not as fast as it could be (the bottleneck being
    ## subsetting an Rle object which needs to be revisited at some point).
    uex_cvg <- revElements(uex_cvg, strand(uex) == "-")

  ## 4) Compute coverage by original exon ('ex_cvg').

    ex2uex <- (seq_along(sm) - cumsum(!is_unique))[sm]  # reverse index
    #stopifnot(identical(ex2uex[uex2ex], seq_along(uex2ex)))  # sanity
    #stopifnot(identical(ex2uex[sm], ex2uex))  # sanity
    #stopifnot(all(uex[ex2uex] == ex))  # sanity

    ex_cvg <- uex_cvg[ex2uex]  # parallel go 'ex'

  ## 5) Compute coverage of each transcript by concatenating coverage of its
  ##    exons.

    ans <- regroup(ex_cvg, transcripts)

  ## 6) Propagate 'mcols(transcripts)'.

    mcols(ans) <- mcols(transcripts)
    ans
}

