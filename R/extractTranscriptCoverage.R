### =========================================================================
### extractTranscriptCoverage()
### -------------------------------------------------------------------------


### A very efficient way to concatenate groups of successive list elements
### in 'x'. 'x' must be a list-like object. 'breakpoints' must be a (possibly
### named) integer vector with no NAs, containing valid indices into 'x' in
### ascending order. If 'breakpoints' has length 0 then 'x' must also have
### length 0. Otherwise 'breakpoints[[length(breakpoints)]]' must be
### 'length(x)'. None of this is checked!
### Return a list-like object of the same elementType() as 'x' but with the
### length and names of 'breakpoints'. Behaves as an endomorphism on a
### CompressedList object.
### Note that
###     .regroup(x, length(x))[[1L]]
### is equivalent to
###     unlist(x, use.names=FALSE)
### TODO: We have this same helper in GenomicFeatures/R/transcripts.R and in
### SplicingGraphs/R/utils.R. Clean this!
.regroup <- function(x, breakpoints)
{
   ans_breakpoints <- end(PartitioningByEnd(x))[breakpoints]
   ans_skeleton <- PartitioningByEnd(ans_breakpoints,
                                     names=names(breakpoints))
   relist(unlist(x, use.names=FALSE), ans_skeleton)
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

    transcripts_breakpoints <- end(PartitioningByEnd(transcripts))
    names(transcripts_breakpoints) <- names(transcripts)
    ans <- .regroup(ex_cvg, transcripts_breakpoints)

  ## 6) Propagate 'mcols(transcripts)'.

    mcols(ans) <- mcols(transcripts)
    ans
}

