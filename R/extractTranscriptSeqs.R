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

#' Extract transcript (or CDS) sequences from chromosome sequences
#'
#' \code{extractTranscriptSeqs} extracts transcript (or CDS) sequences from an
#' object representing a single chromosome or a collection of chromosomes.
#'
#'
#' @aliases extractTranscriptSeqs extractTranscriptSeqs,DNAString-method
#' extractTranscriptSeqs,ANY-method
#' @param x An object representing a single chromosome or a collection of
#' chromosomes.  More precisely, \code{x} can be a \link[Biostrings]{DNAString}
#' object (single chromosome), or a \link[BSgenome]{BSgenome} object
#' (collection of chromosomes).
#'
#' Other objects representing a collection of chromosomes are supported (e.g.
#' \link[Rsamtools]{FaFile} objects in the \pkg{Rsamtools} package) as long as
#' \code{\link[GenomeInfoDb]{seqinfo}} and \code{\link[Biostrings]{getSeq}}
#' work on them.
#' @param transcripts An object representing the exon ranges of each transcript
#' to extract.
#'
#' More precisely: \itemize{ \item If \code{x} is a
#' \link[Biostrings]{DNAString} object, then \code{transcripts} must be an
#' \link[IRanges]{IntegerRangesList} object.
#'
#' \item If \code{x} is a \link[BSgenome]{BSgenome} object or any object
#' representing a collection of chromosomes, then \code{transcripts} must be a
#' \link[GenomicRanges]{GRangesList} object or any object for which
#' \code{\link{exonsBy}} is implemented (e.g. a \link{TxDb} or
#' \link[ensembldb]{EnsDb} object). If the latter, then it's first turned into
#' a \link[GenomicRanges]{GRangesList} object with
#' \code{\link{exonsBy}(transcripts, by="tx", ...)}.  }
#'
#' Note that, for each transcript, the exons must be ordered by ascending
#' \emph{rank}, that is, by ascending position \emph{in the transcript} (when
#' going in the 5' to 3' direction). This generally means (but not always) that
#' they are also ordered from 5' to 3' on the reference genome.  More
#' precisely: \itemize{ \item For a transcript located on the plus strand, the
#' exons will typically (but not necessarily) be ordered by ascending position
#' on the reference genome.  \item For a transcript located on the minus
#' strand, the exons will typically (but not necessarily) be ordered by
#' descending position on the reference genome.  } If \code{transcripts} was
#' obtained with \code{\link{exonsBy}} (see above), then the exons are
#' guaranteed to be ordered by ascending rank. See \code{?\link{exonsBy}} for
#' more information.
#' @param ...  Additional arguments, for use in specific methods.
#'
#' For the default method, additional arguments are allowed only when
#' \code{transcripts} is not a \link[GenomicRanges]{GRangesList} object, in
#' which case they are passed to the internal call to \code{\link{exonsBy}}
#' (see above).
#' @param strand Only supported when \code{x} is a \link[Biostrings]{DNAString}
#' object.
#'
#' Can be an atomic vector, a factor, or an \link[S4Vectors]{Rle} object, in
#' which case it indicates the strand of each transcript (i.e. all the exons in
#' a transcript are considered to be on the same strand).  More precisely: it's
#' turned into a factor (or factor-\link[S4Vectors]{Rle}) that has the
#' "standard strand levels" (this is done by calling the
#' \code{\link[BiocGenerics]{strand}} function on it). Then it's recycled to
#' the length of \link[IRanges]{IntegerRangesList} object \code{transcripts} if
#' needed. In the resulting object, the i-th element is interpreted as the
#' strand of all the exons in the i-th transcript.
#'
#' \code{strand} can also be a list-like object, in which case it indicates the
#' strand of each exon, individually. Thus it must have the same \emph{shape}
#' as \link[IRanges]{IntegerRangesList} object \code{transcripts} (i.e. same
#' length plus \code{strand[[i]]} must have the same length as
#' \code{transcripts[[i]]} for all \code{i}).
#'
#' \code{strand} can only contain \code{"+"} and/or \code{"-"} values.
#' \code{"*"} is not allowed.
#' @return A \link[Biostrings]{DNAStringSet} object \emph{parallel} to
#' \code{transcripts}, that is, the i-th element in it is the sequence of the
#' i-th transcript in \code{transcripts}.
#' @author Hervé Pagès
#' @seealso \itemize{ \item \code{\link{coverageByTranscript}} for computing
#' coverage by transcript (or CDS) of a set of ranges.
#'
#' \item \code{\link{transcriptLengths}} for extracting the transcript lengths
#' (and other metrics) from a \link{TxDb} object.
#'
#' \item \code{\link{extendExonsIntoIntrons}} for extending exons into their
#' adjacent introns.
#'
#' \item The \code{\link{transcriptLocs2refLocs}} function for converting
#' transcript-based locations into reference-based locations.
#'
#' \item The \code{\link[BSgenome]{available.genomes}} function in the
#' \pkg{BSgenome} package for checking avaibility of BSgenome data packages
#' (and installing the desired one).
#'
#' \item The \link[Biostrings]{DNAString} and \link[Biostrings]{DNAStringSet}
#' classes defined and documented in the \pkg{Biostrings} package.
#'
#' \item The \code{\link[Biostrings]{translate}} function in the
#' \pkg{Biostrings} package for translating DNA or RNA sequences into amino
#' acid sequences.
#'
#' \item The \link[GenomicRanges]{GRangesList} class defined and documented in
#' the \pkg{GenomicRanges} package.
#'
#' \item The \link[IRanges]{IntegerRangesList} class defined and documented in
#' the \pkg{IRanges} package.
#'
#' \item The \code{\link{exonsBy}} function for extracting exon ranges grouped
#' by transcript.
#'
#' \item The \link{TxDb} class.  }
#' @keywords manip
#' @examples
#'
#' ## ---------------------------------------------------------------------
#' ## 1. A TOY EXAMPLE
#' ## ---------------------------------------------------------------------
#'
#' library(Biostrings)
#'
#' ## A chromosome of length 30:
#' x <- DNAString("ATTTAGGACACTCCCTGAGGACAAGACCCC")
#'
#' ## 2 transcripts on 'x':
#' tx1 <- IRanges(1, 8)            # 1 exon
#' tx2 <- c(tx1, IRanges(12, 30))  # 2 exons
#' transcripts <- IRangesList(tx1=tx1, tx2=tx2)
#' extractTranscriptSeqs(x, transcripts)
#'
#' ## By default, all the exons are considered to be on the plus strand.
#' ## We can use the 'strand' argument to tell extractTranscriptSeqs()
#' ## to extract them from the minus strand.
#'
#' ## Extract all the exons from the minus strand:
#' extractTranscriptSeqs(x, transcripts, strand="-")
#'
#' ## Note that, for a transcript located on the minus strand, the exons
#' ## should typically be ordered by descending position on the reference
#' ## genome in order to reflect their rank in the transcript:
#' extractTranscriptSeqs(x, IRangesList(tx1=tx1, tx2=rev(tx2)), strand="-")
#'
#' ## Extract the exon of the 1st transcript from the minus strand:
#' extractTranscriptSeqs(x, transcripts, strand=c("-", "+"))
#'
#' ## Extract the 2nd exon of the 2nd transcript from the minus strand:
#' extractTranscriptSeqs(x, transcripts, strand=list("-", c("+", "-")))
#'
#' ## ---------------------------------------------------------------------
#' ## 2. A REAL EXAMPLE
#' ## ---------------------------------------------------------------------
#'
#' ## Load a genome:
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome <- BSgenome.Hsapiens.UCSC.hg19
#'
#' ## Load a TxDb object:
#' txdb_file <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                          package="GenomicFeatures")
#' txdb <- loadDb(txdb_file)
#'
#' ## Check that 'txdb' is based on the hg19 assembly:
#' txdb
#'
#' ## Extract the exon ranges grouped by transcript from 'txdb':
#' transcripts <- exonsBy(txdb, by="tx", use.names=TRUE)
#'
#' ## Extract the transcript sequences from the genome:
#' tx_seqs <- extractTranscriptSeqs(genome, transcripts)
#' tx_seqs
#'
#' ## A sanity check:
#' stopifnot(identical(width(tx_seqs), unname(sum(width(transcripts)))))
#'
#' ## Note that 'tx_seqs' can also be obtained with:
#' extractTranscriptSeqs(genome, txdb, use.names=TRUE)
#'
#' ## ---------------------------------------------------------------------
#' ## 3. USING extractTranscriptSeqs() TO EXTRACT CDS SEQUENCES
#' ## ---------------------------------------------------------------------
#'
#' cds <- cdsBy(txdb, by="tx", use.names=TRUE)
#' cds_seqs <- extractTranscriptSeqs(genome, cds)
#' cds_seqs
#'
#' ## A sanity check:
#' stopifnot(identical(width(cds_seqs), unname(sum(width(cds)))))
#'
#' ## Note that, alternatively, the CDS sequences can be obtained from the
#' ## transcript sequences by removing the 5' and 3' UTRs:
#' tx_lens <- transcriptLengths(txdb, with.utr5_len=TRUE, with.utr3_len=TRUE)
#' stopifnot(identical(tx_lens$tx_name, names(tx_seqs)))  # sanity
#' ## Keep the rows in 'tx_lens' that correspond to a sequence in 'cds_seqs'
#' ## and put them in the same order as in 'cds_seqs':
#' m <- match(names(cds_seqs), names(tx_seqs))
#' tx_lens <- tx_lens[m, ]
#' utr5_width <- tx_lens$utr5_len
#' utr3_width <- tx_lens$utr3_len
#' cds_seqs2 <- narrow(tx_seqs[m],
#'                     start=utr5_width+1L, end=-(utr3_width+1L))
#' stopifnot(identical(as.character(cds_seqs2), as.character(cds_seqs)))
#'
#' ## ---------------------------------------------------------------------
#' ## 4. TRANSLATE THE CDS SEQUENCES
#' ## ---------------------------------------------------------------------
#'
#' prot_seqs <- translate(cds_seqs, if.fuzzy.codon="solve")
#'
#' ## Note that, by default, translate() uses The Standard Genetic Code to
#' ## translate codons into amino acids. However, depending on the organism,
#' ## a different genetic code might be needed to translate CDS sequences
#' ## located on the mitochodrial chromosome. For example, for vertebrates,
#' ## the following code could be used to correct 'prot_seqs':
#' SGC1 <- getGeneticCode("SGC1")
#' chrM_idx <- which(all(seqnames(cds) == "chrM"))
#' prot_seqs[chrM_idx] <- translate(cds_seqs[chrM_idx], genetic.code=SGC1,
#'                                  if.fuzzy.codon="solve")
#'
#' @export
setGeneric("extractTranscriptSeqs", signature="x",
    function(x, transcripts, ...) standardGeneric("extractTranscriptSeqs")
)

setMethod("extractTranscriptSeqs", "DNAString",
    function(x, transcripts, strand="+")
    {
        if (!is(transcripts, "IntegerRangesList"))
            stop(wmsg("when 'x' is a DNAString object, ",
                      "'transcripts' must be an IntegerRangesList object"))
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
    target <- sequence(tx1_eltNROWS, from=min_rank)
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

.extractTranscriptSeqs_default <- function(x, transcripts, ...)
{
    if (is(transcripts, "GRangesList")) {
        if (length(list(...)) != 0L)
            stop(wmsg("additional arguments are allowed only when ",
                      "'transcripts' is not a GRangesList object"))
    } else {
        transcripts <- try(exonsBy(transcripts, by="tx", ...),
                           silent=TRUE)
        if (is(transcripts, "try-error"))
            stop(wmsg("failed to extract the exon ranges from 'transcripts' ",
                      "with exonsBy(transcripts, by=\"tx\", ...)"))
    }
    idx1 <- which(elementNROWS(transcripts) != 0L)
    tx1 <- transcripts[idx1]
    .check_exon_chrom(tx1)
    .check_exon_rank(tx1)

    tx1_seqlevels_in_use <- seqlevelsInUse(tx1)
    x_seqlevels <- seqlevels(x)
    ok <- tx1_seqlevels_in_use %in% x_seqlevels
    if (!all(ok)) {
        if (all(!ok))
            stop(wmsg("the transcripts in 'transcripts' are on chromosomes ",
                      "that are not in 'x'"))
        seqlevel_not_in_x <- tx1_seqlevels_in_use[!ok][[1L]]
        stop(wmsg("some transcripts in 'transcripts' are on chromosomes ",
                  "that are not in 'x' (e.g. some transcripts are on ",
                  "chromosome \"", seqlevel_not_in_x, "\" but this ",
                  "chromosome is not in 'x')"))
    }
    seqlevels(tx1) <- tx1_seqlevels_in_use
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
setMethod("extractTranscriptSeqs", "ANY", .extractTranscriptSeqs_default)

