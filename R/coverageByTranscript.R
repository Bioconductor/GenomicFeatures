### =========================================================================
### coverageByTranscript()
### -------------------------------------------------------------------------


### Infers the seqlengths based on the ranges located on each sequence.
### 'x' can be a GRanges, GRangesList, GAlignments, GAlignmentPairs, or
### GAlignmentsList object. Maybe .infer_seqlengths() should be a generic
### function?
### Return an integer vector parallel to 'seqlevels(x)' with NAs where the
### sequence length could not be inferred because the sequence has no range
### on it.
.infer_seqlengths <- function(x)
{
    if (is(x, "CompressedList") || is(x, "GAlignmentPairs"))
        x <- unlist(x, use.names=FALSE)
    ## Now 'x' should only be a GRanges or GAlignments object.
    pmax(as.integer(end(range(split(ranges(x), seqnames(x))))), 0L)
}

### Infer seqlengths that are missing from 'x' and 'transcripts'.
.merge_seqinfo_and_infer_missing_seqlengths <- function(x, transcripts)
{
    ans <- merge(seqinfo(x), seqinfo(transcripts))
    seqlevels(x) <- seqlevels(transcripts) <- seqlevels(ans)
    x_inferred_seqlengths <- .infer_seqlengths(x)
    transcripts_inferred_seqlengths <- .infer_seqlengths(transcripts)
    ## 'x_inferred_seqlengths' and 'transcripts_inferred_seqlengths' are
    ## guaranteed to be parallel.
    inferred_seqlengths <- pmax(x_inferred_seqlengths,
                                transcripts_inferred_seqlengths,
                                na.rm=TRUE)
    missing_idx <- which(is.na(seqlengths(ans)))
    seqlengths(ans)[missing_idx] <- inferred_seqlengths[missing_idx]
    ans
}

### Computing coverage by transcript (or CDS) of a set of ranges is an
### operation that feels a lot like extracting transcript sequences from a
### genome. Defined as an ordinary function for now.


#' Compute coverage by transcript (or CDS) of a set of ranges
#' 
#' \code{coverageByTranscript} computes the transcript (or CDS) coverage of a
#' set of ranges.
#' 
#' \code{pcoverageByTranscript} is a version of \code{coverageByTranscript}
#' that operates element-wise.
#' 
#' 
#' @aliases coverageByTranscript pcoverageByTranscript
#' @param x An object representing a set of ranges (typically aligned reads).
#' \link[GenomicRanges]{GRanges}, \link[GenomicRanges]{GRangesList},
#' \link[GenomicAlignments]{GAlignments},
#' \link[GenomicAlignments]{GAlignmentPairs}, and
#' \link[GenomicAlignments]{GAlignmentsList} objects are supported.
#' 
#' More generally, for \code{coverageByTranscript} \code{x} can be any object
#' for which \code{\link[GenomeInfoDb]{seqinfo}()} and
#' \code{\link[GenomicRanges]{coverage}()} are supported (e.g. a
#' \link[Rsamtools]{BamFile} object).  Note that, for such objects,
#' \code{coverage()} is expected to return an \link[IRanges]{RleList} object
#' whose names are \code{seqlevels(x)}).
#' 
#' More generally, for \code{pcoverageByTranscript} \code{x} can be any object
#' for which \code{\link[GenomicRanges]{grglist}()} is supported.  It should
#' have the length of \code{transcripts} or length 1. If the latter, it is
#' recycled to the length of \code{transcripts}.
#' @param transcripts A \link[GenomicRanges]{GRangesList} object representing
#' the exons of each transcript for which to compute coverage. For each
#' transcript, the exons must be ordered by \emph{ascending rank}, that is, by
#' their position in the transcript. This means that, for a transcript located
#' on the minus strand, the exons should typically be ordered by descending
#' position on the reference genome. If \code{transcripts} was obtained with
#' \code{\link{exonsBy}}, then the exons are guaranteed to be ordered by
#' ascending rank. See \code{?\link{exonsBy}} for more information.
#' 
#' Alternatively, \code{transcripts} can be a \link{TxDb} object, or any
#' \link{TxDb}-like object that supports the \code{\link{exonsBy}()} extractor
#' (e.g. an \link[ensembldb]{EnsDb} object). In this case it is replaced with
#' the \link[GenomicRanges]{GRangesList} object returned by
#' \code{\link{exonsBy}(transcripts, by="tx", use.names=TRUE)}.
#' 
#' For \code{pcoverageByTranscript}, \code{transcripts} should have the length
#' of \code{x} or length 1. If the latter, it is recycled to the length of
#' \code{x}.
#' @param ignore.strand TRUE or FALSE. If FALSE (the default) then the strand
#' of a range in \code{x} and exon in \code{transcripts} must be the same in
#' order for the range to contribute coverage to the exon. If TRUE then the
#' strand is ignored.
#' @param ...  Additional arguments passed to the internal call to
#' \code{\link[GenomicRanges]{grglist}()}.  More precisely, when \code{x} is
#' not a \link[GenomicRanges]{GRanges} or \link[GenomicRanges]{GRangesList}
#' object, \code{pcoverageByTranscript} replace it with the
#' \link[GenomicRanges]{GRangesList} object returned by
#' \code{\link[GenomicRanges]{grglist}(x, ...)}.
#' @return An \link[IRanges]{RleList} object \emph{parallel} to
#' \code{transcripts}, that is, the i-th element in it is an
#' integer-\link[S4Vectors]{Rle} representing the coverage of the i-th
#' transcript in \code{transcripts}.  Its \code{lengths()} is guaranteed to be
#' identical to \code{sum(width(transcripts))}. The names and metadata columns
#' on \code{transcripts} are propagated to it.
#' @author Hervé Pagès
#' @seealso \itemize{ \item \code{\link{transcripts}},
#' \code{\link{transcriptsBy}}, and \code{\link{transcriptsByOverlaps}}, for
#' extracting genomic feature locations from a \link{TxDb}-like object.
#' 
#' \item \code{\link{transcriptLengths}} for extracting the transcript lengths
#' (and other metrics) from a \link{TxDb} object.
#' 
#' \item \code{\link{extractTranscriptSeqs}} for extracting transcript (or CDS)
#' sequences from chromosome sequences.
#' 
#' \item The \link[IRanges]{RleList} class defined and documented in the
#' \pkg{IRanges} package.
#' 
#' \item The \link[GenomicRanges]{GRangesList} class defined and documented in
#' the \pkg{GenomicRanges} package.
#' 
#' \item The \code{\link[GenomicRanges]{coverage}} methods defined in the
#' \pkg{GenomicRanges} package.
#' 
#' \item The \code{\link{exonsBy}} function for extracting exon ranges grouped
#' by transcript.
#' 
#' \item \code{\link[GenomicAlignments]{findCompatibleOverlaps}} in the
#' \pkg{GenomicAlignments} package for finding which reads are
#' \emph{compatible} with the splicing of which transcript.  }
#' @keywords manip
#' @examples
#' 
#' ## ---------------------------------------------------------------------
#' ## 1. A SIMPLE ARTIFICIAL EXAMPLE WITH ONLY ONE TRANSCRIPT
#' ## ---------------------------------------------------------------------
#' 
#' ## Get some transcripts:
#' library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
#' txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
#' dm3_transcripts <- exonsBy(txdb, by="tx", use.names=TRUE)
#' dm3_transcripts
#' 
#' ## Let's pick up the 1st transcript: FBtr0300689. It as 2 exons and 1
#' ## intron:
#' my_transcript <- dm3_transcripts["FBtr0300689"]
#' 
#' ## Let's create 3 artificial aligned reads. We represent them as a
#' ## GRanges object of length 3 that contains the genomic positions of
#' ## the 3 reads. Note that these reads are simple alignments i.e. each
#' ## of them can be represented with a single range. This would not be
#' ## the case if they were junction reads.
#' my_reads <- GRanges(c("chr2L:7531-7630",
#'                       "chr2L:8101-8200",
#'                       "chr2L:8141-8240"))
#' 
#' ## The coverage of the 3 reads on the reference genome is:
#' coverage(my_reads)
#' 
#' ## As you can see, all the genomic positions in the 3 ranges participate
#' ## to the coverage. This can be confirmed by comparing:
#' sum(coverage(my_reads))
#' ## with:
#' sum(width(my_reads))
#' ## They should always be the same.
#' 
#' ## When computing the coverage on a transcript, only the part of the
#' ## read that overlaps with the transcript participates to the coverage.
#' ## Let's look at the individual coverage of each read on transcript
#' ## FBtr0300689:
#' 
#' ## The 1st read is fully contained within the 1st exon:
#' coverageByTranscript(my_reads[1], my_transcript)
#' 
#' ## Note that the length of the Rle (1880) is the length of the transcript.
#' 
#' ## The 2nd and 3rd reads overlap the 2 exons and the intron. Only the
#' ## parts that overlap the exons participate to coverage:
#' coverageByTranscript(my_reads[2], my_transcript)
#' coverageByTranscript(my_reads[3], my_transcript)
#' 
#' ## The coverage of the 3 reads together is:
#' coverageByTranscript(my_reads, my_transcript)
#' 
#' ## Note that this is the sum of the individual coverages. This can be
#' ## checked with:
#' stopifnot(all(
#'   coverageByTranscript(my_reads, my_transcript)
#'   ==
#'   Reduce("+", lapply(seq_along(my_reads),
#'       function(i) coverageByTranscript(my_reads[i], my_transcript)), 0L)
#' ))
#' 
#' ## ---------------------------------------------------------------------
#' ## 2. COMPUTE THE FULL TRANSCRIPTOME COVERAGE OF A SET OF ALIGNED READS
#' ## ---------------------------------------------------------------------
#' 
#' ## Load the aligned reads:
#' library(pasillaBamSubset)
#' library(GenomicAlignments)
#' reads <- readGAlignments(untreated1_chr4())
#' 
#' ## Compute the full transcriptome coverage by calling
#' ## coverageByTranscript() on 'dm3_transcripts':
#' tx_cvg <- coverageByTranscript(reads, dm3_transcripts, ignore.strand=TRUE)
#' tx_cvg
#' 
#' ## A sanity check:
#' stopifnot(identical(lengths(tx_cvg), sum(width(dm3_transcripts))))
#' 
#' ## We can also use pcoverageByTranscript() to compute 'tx_cvg'.
#' ## For this we first create a GAlignmentsList object "parallel" to
#' ## 'dm3_transcripts' where the i-th list element contains the aligned
#' ## reads that overlap with the i-th transcript:
#' hits <- findOverlaps(reads, dm3_transcripts, ignore.strand=TRUE)
#' tx2reads <- setNames(as(t(hits), "List"), names(dm3_transcripts))
#' reads_by_tx <- extractList(reads, tx2reads)  # GAlignmentsList object
#' reads_by_tx
#' 
#' ## Call pcoverageByTranscript():
#' tx_cvg2 <- pcoverageByTranscript(reads_by_tx, dm3_transcripts,
#'                                  ignore.strand=TRUE)
#' stopifnot(identical(tx_cvg, tx_cvg2))
#' 
#' ## A more meaningful coverage is obtained by counting for each
#' ## transcript only the reads that are *compatible* with its splicing:
#' compat_hits <- findCompatibleOverlaps(reads, dm3_transcripts)
#' tx2reads <- setNames(as(t(compat_hits), "List"), names(dm3_transcripts))
#' compat_reads_by_tx <- extractList(reads, tx2reads)
#' 
#' tx_compat_cvg <- pcoverageByTranscript(compat_reads_by_tx,
#'                                        dm3_transcripts,
#'                                        ignore.strand=TRUE)
#' ## A sanity check:
#' stopifnot(all(all(tx_compat_cvg <= tx_cvg)))
#' 
#' ## ---------------------------------------------------------------------
#' ## 3. COMPUTE CDS COVERAGE OF A SET OF ALIGNED READS
#' ## ---------------------------------------------------------------------
#' 
#' ## coverageByTranscript() can also be used to compute CDS coverage:
#' cds <- cdsBy(txdb, by="tx", use.names=TRUE)
#' cds_cvg <- coverageByTranscript(reads, cds, ignore.strand=TRUE)
#' cds_cvg
#' 
#' ## A sanity check:
#' stopifnot(identical(lengths(cds_cvg), sum(width(cds))))
#' 
#' ## ---------------------------------------------------------------------
#' ## 4. ALTERNATIVELY, THE CDS COVERAGE CAN BE OBTAINED FROM THE
#' ##    TRANSCRIPT COVERAGE BY TRIMMING THE 5' AND 3' UTRS
#' ## ---------------------------------------------------------------------
#' 
#' tx_lens <- transcriptLengths(txdb, with.utr5_len=TRUE, with.utr3_len=TRUE)
#' stopifnot(identical(tx_lens$tx_name, names(tx_cvg)))  # sanity
#' 
#' ## Keep the rows in 'tx_lens' that correspond to a list element in
#' ## 'cds_cvg' and put them in the same order as in 'cds_cvg':
#' m <- match(names(cds_cvg), names(tx_cvg))
#' tx_lens <- tx_lens[m, ]
#' utr5_width <- tx_lens$utr5_len
#' utr3_width <- tx_lens$utr3_len
#' cds_cvg2 <- windows(tx_cvg[m], start=1L+utr5_width, end=-1L-utr3_width)
#' 
#' ## A sanity check:
#' stopifnot(identical(cds_cvg2, cds_cvg))
#' 
#' @export coverageByTranscript
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
    if (!isTRUEorFALSE(ignore.strand))
        stop(wmsg("'ignore.strand' must be TRUE or FALSE"))

  ## STEP 1 - Infer seqlengths that are missing from 'x' and 'transcripts'.

    ## This will guarantee that subsetting the named RleList object
    ## representing the read coverage by the GRanges object representing
    ## the set of unique exons won't fail due to out-of-bounds ranges in
    ## the subscript. See STEP 3 below.
    seqinfo(x) <- .merge_seqinfo_and_infer_missing_seqlengths(x, transcripts)

  ## STEP 2 - Compute unique exons ('uex').

    ex <- unlist(transcripts, use.names=FALSE)
    ## We could simply do 'uex <- unique(ex)' here but we're going to need
    ## 'sm' and 'is_unique' later to compute the "reverse index" so we compute
    ## them now and use them to extract the unique exons. That way we hash
    ## 'ex' only once (the expensive operation).
    sm <- selfmatch(ex)  # uses a hash table internally
    is_unique <- sm == seq_along(sm)
    uex2ex <- which(is_unique)  # index of unique exons
    uex <- ex[uex2ex]  # unique exons

  ## STEP 3 - Compute coverage for each unique exon ('uex_cvg').

    #There doesn't seem to be much benefit in doing this.
    #x <- subsetByOverlaps(x, transcripts, ignore.strand=TRUE)
    if (ignore.strand) {
        cvg <- coverage(x)
        ## Because we've inferred the seqlengths that are missing from 'x'
        ## and 'transcripts' (see STEP 1 above), 'uex' should never contain
        ## out-of-bounds ranges i.e. the list elements in 'cvg' should always
        ## be Rle objects that are long enough with respect to the ranges
        ## in 'uex'.
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
        ## Because we've inferred the seqlengths that are missing from 'x'
        ## and 'transcripts' (see STEP 1 above), 'uex' should never contain
        ## out-of-bounds ranges i.e. the list elements in 'cvg1' and 'cvg2'
        ## should always be Rle objects that are long enough with respect to
        ## the ranges in 'uex'.
        uex_cvg <- cvg1[uex]
        uex_cvg[is_minus_ex] <- cvg2[uex[is_minus_ex]]
    }

  ## STEP 4 - Flip coverage for exons on minus strand.

    ## It feels like this is not as fast as it could be (the bottleneck being
    ## subsetting an Rle object which needs to be revisited at some point).
    uex_cvg <- revElements(uex_cvg, strand(uex) == "-")

  ## STEP 5 - Compute coverage by original exon ('ex_cvg').

    ex2uex <- (seq_along(sm) - cumsum(!is_unique))[sm]  # reverse index
    #stopifnot(identical(ex2uex[uex2ex], seq_along(uex2ex)))  # sanity
    #stopifnot(identical(ex2uex[sm], ex2uex))  # sanity
    #stopifnot(all(uex[ex2uex] == ex))  # sanity

    ex_cvg <- uex_cvg[ex2uex]  # parallel go 'ex'

  ## STEP 6 - Compute coverage of each transcript by concatenating coverage of
  ##          its exons.

    ans <- IRanges:::regroupBySupergroup(ex_cvg, transcripts)

  ## STEP 7 - Propagate 'mcols(transcripts)'.

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

