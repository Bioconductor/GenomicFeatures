### =========================================================================
### coverageByTranscript()
### -------------------------------------------------------------------------


### Computing coverage by transcript (or CDS) of a set of ranges is an
### operation that feels a lot like extracting transcript sequences from a
### genome. Defined as an ordinary function for now.
coverageByTranscript <- function(x, transcripts, ignore.strand=FALSE)
{
    if (!is(transcripts, "GRangesList")) {
        transcripts <- try(exonsBy(transcripts, by="tx", use.names=TRUE),
                           silent=TRUE)
        if (is(transcripts, "try-error"))
            stop(wmsg("failed to extract the exon ranges ",
                      "from 'transcripts' with ",
                      "exonsBy(transcripts, by=\"tx\", use.names=TRUE)"))
    }
    seqinfo(x) <- merge(seqinfo(x), seqinfo(transcripts))
    if (!isTRUEorFALSE(ignore.strand))
        stop(wmsg("'ignore.strand' must be TRUE or FALSE"))

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
    #x <- subsetByOverlaps(x, transcripts, ignore.strand=TRUE)
    if (ignore.strand) {
        cvg <- coverage(x)
        uex_cvg <- cvg[uex]  # parallel to 'uex'
    } else {
        x1 <- x[strand(x) %in% c("+", "*")]
        x2 <- x[strand(x) %in% c("-", "*")]
        cvg1 <- coverage(x1)
        cvg2 <- coverage(x2)
        is_plus_ex <- strand(uex) == "+"
        is_minus_ex <- strand(uex) == "-"
        if (!identical(is_plus_ex, !is_minus_ex))
            stop(wmsg("'transcripts' has exons on the * strand. ",
                      "This is not supported at the moment."))
        uex_cvg <- cvg1[uex]
        uex_cvg[is_minus_ex] <- cvg2[uex[is_minus_ex]]
    }

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

    ans <- IRanges:::regroupBySupergroup(ex_cvg, transcripts)

  ## 6) Propagate 'mcols(transcripts)'.

    mcols(ans) <- mcols(transcripts)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### projectOnTranscripts()
###
### NOT exported! Used in pcoverageByTranscript() below.
###

### Transform genome-based coordinates into transcriptome-based coordinates.
### 'transcripts' must be a GRanges object representing the exons of a single
### transcript or a GRangesList object of exons grouped by transcript. If the
### former then 'start' and 'end' must be numeric vectors parallel to
### 'transcripts' and the function returns 2 numeric vectors also parallel to
### 'transcripts'. If the latter then they must be NumericList objects with
### the shape of 'transcripts' and the function returns 2 NumericList objects
### also with the shape of 'transcripts'.
.genome2txcoords <- function(transcripts, start, end=start)
{
    offset_in_exon <- start - start(transcripts)
    offset_in_exon2 <- end(transcripts) - end
    cumsumwtx <- cumsum(width(transcripts))
    end_in_tx <- cumsumwtx - offset_in_exon2
    end_in_tx2 <- cumsumwtx - offset_in_exon
    strand_is_minus <- unname(strand(transcripts) == "-")
    if (is(transcripts, "GRangesList")) {
        strand_is_minus <- unlist(strand_is_minus, use.names=FALSE)
        tmp <- extractROWS(unlist(offset_in_exon2, use.names=FALSE),
                           strand_is_minus)
        offset_in_exon <- relist(
               replaceROWS(unlist(offset_in_exon, use.names=FALSE),
                           strand_is_minus, tmp),
                           offset_in_exon)
        tmp <- extractROWS(unlist(end_in_tx2, use.names=FALSE),
                           strand_is_minus)
        end_in_tx <- relist(
               replaceROWS(unlist(end_in_tx, use.names=FALSE),
                           strand_is_minus, tmp),
                           end_in_tx)
    } else {  # GRanges
        tmp <- extractROWS(offset_in_exon2, strand_is_minus)
        offset_in_exon <- replaceROWS(offset_in_exon, strand_is_minus, tmp)
        tmp <- extractROWS(end_in_tx2, strand_is_minus)
        end_in_tx <- replaceROWS(end_in_tx, strand_is_minus, tmp)
    }
    list(end_in_tx=end_in_tx, offset_in_exon=offset_in_exon)
}

### Put nice transcriptome costume on.
.set_transcriptome_seqinfo <- function(x, transcripts)
{
    stopifnot(is(x, "GRangesList"), is(transcripts, "GRangesList"))
    stopifnot(length(transcripts) == length(x) || length(transcripts) == 1L)

    ## Temporarily put all the ranges in 'x' on the 1st level and drop
    ## all other levels.
    seqlevel1 <- seqlevels(x)[[1L]]
    tmp_seqnames <- relist(Rle(factor(seqlevel1, levels=seqlevels(x)),
                               length(unlist(x, use.names=FALSE))), x)
    suppressWarnings(seqnames(x) <- tmp_seqnames)
    seqlevels(x) <- seqlevel1

    ## Prepare transcriptome-based seqlevels.
    ans_seqlevels <- names(transcripts)
    if (is.null(ans_seqlevels)) {
        seqlevels_are_bad <- TRUE
    } else {
        seqlevels_are_bad <- anyDuplicated(ans_seqlevels) ||
                             any(ans_seqlevels %in% c(NA_character_, ""))
        if (seqlevels_are_bad)
            warning(wmsg("The names on 'transcripts' could not be used to set ",
                         "the seqlevels of the returned object (because they ",
                         "contain duplicates, NAs, and/or empty strings). ",
                         "Setting artificial seqlevels TX1, TX2, etc..."))
    }
    if (seqlevels_are_bad)
        ans_seqlevels <- paste0("TX", as.character(seq_along(transcripts)))

    ## Set transcriptome-based seqinfo.
    ans_seqinfo <- Seqinfo(ans_seqlevels,
                           unname(sum(width(transcripts))),
                           logical(length(transcripts)),
                           genome(transcripts))
    new2old <- rep.int(NA_integer_, length(ans_seqinfo))
    new2old[1L] <- 1L
    suppressWarnings(seqinfo(x, new2old=new2old) <- ans_seqinfo)

    ## Set transcriptome-based seqnames.
    if (length(ans_seqlevels) != 1L) {
        ans_seqnames <- relist(Rle(factor(ans_seqlevels, levels=ans_seqlevels),
                                   elementNROWS(x)), x)
        seqnames(x) <- ans_seqnames
    }

    ## Set strand to "*".
    strand(x) <- "*"
    x
}

### If 'keep.nohit.exons' is TRUE, return a GRangesList object with the same
### shape as 'transcripts'. Names and metadata columns on 'transcripts' are
### propagated.
.project_GRanges_on_transcripts <- function(x, transcripts,
                                            ignore.strand=FALSE,
                                            keep.nohit.exons=FALSE,
                                            set.transcriptome.seqinfo=TRUE)
{
    stopifnot(is(x, "GRanges"), is(transcripts, "GRangesList"))

    ## Recycle arguments.
    transcripts_was_recycled <- FALSE
    if (length(x) != length(transcripts)) {
        if (length(x) == 1L) {
            x <- rep.int(x, length(transcripts))
        } else if (length(transcripts) == 1L) {
            transcripts <- rep.int(transcripts, length(x))
            transcripts_was_recycled <- TRUE
        } else {
            stop(wmsg("when 'x' and 'transcripts' don't have the ",
                      "same length, one of them must have length 1"))
        }
    }
    if (!isTRUEorFALSE(keep.nohit.exons))
        stop(wmsg("'keep.nohit.exons' must be TRUE or FALSE"))

    pint <- pintersect(transcripts, x, ignore.strand=ignore.strand)

    ## Transform genome-based coordinates in 'pint' into transcriptome-based
    ## coordinates.
    txcoords <- .genome2txcoords(transcripts, start(pint), end(pint))
    seqlengths(pint)[] <- NA_integer_  # so no out-of-bound range when shifting
    ans <- shift(pint, txcoords$end_in_tx - end(pint))

    ## Add "offset_in_exon" inner metadata column to 'ans'.
    unlisted_ans <- unlist(ans, use.names=FALSE)
    mcols(unlisted_ans)$offset_in_exon <- unlist(txcoords$offset_in_exon,
                                                 use.names=FALSE)

    ## 'pint' typically contains many "nohit ranges" that are zero-width
    ## ranges artificially created by pintersect(), one for each exon in
    ## 'transcripts' that receives no hit. This allows 'pint' to have the
    ## same shape as 'transcripts' which is required for .genome2txcoords()
    ## to work. If 'keep.nohit.exons' is FALSE then we remove those "nohit
    ## ranges" from 'ans'.
    if (!keep.nohit.exons) {
        is_hit <- relist(mcols(unlisted_ans)$hit, ans)
        mcols(unlisted_ans)$hit <- NULL
        ans <- relist(unlisted_ans, ans)[is_hit]
    } else {
        ans <- relist(unlisted_ans, ans)
    }

    if (!set.transcriptome.seqinfo)
        return(ans)

    ## Put nice transcriptome costume on.
    if (transcripts_was_recycled)
        transcripts <- transcripts[1L]
    .set_transcriptome_seqinfo(ans, transcripts)
}

### Return a GRangesList object parallel to 'transcripts' (but the shape is
### different). Names and metadata columns on 'transcripts' are propagated.
.project_GRangesList_on_transcripts <- function(x, transcripts,
                                                ignore.strand=FALSE,
                                                keep.nohit.exons=FALSE)
{
    stopifnot(is(x, "GRangesList"), is(transcripts, "GRangesList"))

    ## Recycle arguments.
    transcripts_was_recycled <- FALSE
    if (length(x) != length(transcripts)) {
        if (length(x) == 1L) {
            x <- rep.int(x, length(transcripts))
        } else if (length(transcripts) == 1L) {
            transcripts <- rep.int(transcripts, length(x))
            transcripts_was_recycled <- TRUE
        } else {
            stop(wmsg("when 'x' and 'transcripts' don't have the ",
                      "same length, one of them must have length 1"))
        }
    }

    x_eltNROWS <- elementNROWS(x)
    transcripts2 <- rep.int(transcripts, x_eltNROWS)
    unlisted_x <- unlist(x, use.names=FALSE)
    y <- .project_GRanges_on_transcripts(unlisted_x, transcripts2,
                                         ignore.strand, keep.nohit.exons,
                                         set.transcriptome.seqinfo=FALSE)
    ans <- IRanges:::regroupBySupergroup(y, x)

    ## Put nice transcriptome costume on.
    if (transcripts_was_recycled)
        transcripts <- transcripts[1L]
    .set_transcriptome_seqinfo(ans, transcripts)
}

### 'x' and 'transcripts' must have the same length. Perform "parallel
### projection" of the ranges in 'x' on the exons in 'transcripts'.
setGeneric("projectOnTranscripts", signature="x",
    function(x, transcripts, ...) standardGeneric("projectOnTranscripts")
)

setMethod("projectOnTranscripts", "GRanges",
    .project_GRanges_on_transcripts
)

setMethod("projectOnTranscripts", "GRangesList",
    .project_GRangesList_on_transcripts
)

setMethod("projectOnTranscripts", "ANY",
    function(x, transcripts, ignore.strand=FALSE, keep.nohit.exons=FALSE, ...)
    {
        x <- grglist(x, ...)
        callGeneric()
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pcoverageByTranscript()
###
### "Parallel" version of coverageByTranscript(). 'x' and 'transcripts'
### must have the same length. 'x' can be any object supported by
### projectOnTranscripts().
###

pcoverageByTranscript <- function(x, transcripts, ignore.strand=FALSE, ...)
{
    if (!is(transcripts, "GRangesList")) {
        transcripts <- try(exonsBy(transcripts, by="tx", use.names=TRUE),
                           silent=TRUE)
        if (is(transcripts, "try-error"))
            stop(wmsg("failed to extract the exon ranges ",
                      "from 'transcripts' with ",
                      "exonsBy(transcripts, by=\"tx\", use.names=TRUE)"))
    }
    ranges_on_tx <- projectOnTranscripts(x, transcripts,
                                         ignore.strand=ignore.strand,
                                         keep.nohit.exons=FALSE,
                                         ...)
    ans <- as(coverage(ranges_on_tx), "CompressedRleList")
    names(ans) <- names(transcripts)
    mcols(ans) <- mcols(transcripts)
    ans
}

