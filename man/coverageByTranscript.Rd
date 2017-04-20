\name{coverageByTranscript}

\alias{coverageByTranscript}
\alias{pcoverageByTranscript}

\title{Compute coverage by transcript (or CDS) of a set of ranges}

\description{
  \code{coverageByTranscript} computes the transcript (or CDS) coverage
  of a set of ranges.

  \code{pcoverageByTranscript} is a version of \code{coverageByTranscript}
  that operates element-wise.
}

\usage{
coverageByTranscript(x, transcripts, ignore.strand=FALSE)

pcoverageByTranscript(x, transcripts, ignore.strand=FALSE, ...)
}

\arguments{
  \item{x}{
    An object representing a set of ranges (typically aligned reads).
    \link[GenomicRanges]{GRanges}, \link[GenomicRanges]{GRangesList},
    \link[GenomicAlignments]{GAlignments},
    \link[GenomicAlignments]{GAlignmentPairs}, and
    \link[GenomicAlignments]{GAlignmentsList} objects are supported.

    More generally, for \code{coverageByTranscript} \code{x} can be
    any object for which \code{\link[GenomeInfoDb]{seqinfo}()} and
    \code{\link[GenomicRanges]{coverage}()} are supported (e.g. a
    \link[Rsamtools]{BamFile} object).
    Note that, for such objects, \code{coverage()} is expected to return an
    \link[IRanges]{RleList} object whose names are \code{seqlevels(x)}).

    More generally, for \code{pcoverageByTranscript} \code{x} can be
    any object for which \code{\link[GenomicRanges]{grglist}()} is supported.
    It should have the length of \code{transcripts} or length 1. If the
    latter, it is recycled to the length of \code{transcripts}.
  }
  \item{transcripts}{
    A \link[GenomicRanges]{GRangesList} object representing the exons of
    each transcript for which to compute coverage. For each transcript, the
    exons must be ordered by \emph{ascending rank}, that is, by their position
    in the transcript. This means that, for a transcript located on the minus
    strand, the exons should typically be ordered by descending position on
    the reference genome. If \code{transcripts} was obtained with
    \code{\link{exonsBy}}, then the exons are guaranteed to be ordered by
    ascending rank. See \code{?\link{exonsBy}} for more information.

    Alternatively \code{transcripts} can be any object for which
    \code{\link{exonsBy}} is implemented (e.g. a \link{TxDb} object), in
    which case it is replaced with the \link[GenomicRanges]{GRangesList} object
    returned by \code{\link{exonsBy}(transcripts, by="tx", use.names=TRUE)}.

    For \code{pcoverageByTranscript}, \code{transcripts} should have the
    length of \code{x} or length 1. If the latter, it is recycled to the
    length of \code{x}.
  }
  \item{ignore.strand}{
    TRUE or FALSE. If FALSE (the default) then the strand of a range in
    \code{x} and exon in \code{transcripts} must be the same in order for
    the range to contribute coverage to the exon. If TRUE then the strand
    is ignored.
  }
  \item{...}{
    Additional arguments passed to the internal call to
    \code{\link[GenomicRanges]{grglist}()}.
    More precisely, when \code{x} is not a \link[GenomicRanges]{GRanges}
    or \link[GenomicRanges]{GRangesList} object,
    \code{pcoverageByTranscript} replace it with the
    \link[GenomicRanges]{GRangesList} object returned by
    \code{\link[GenomicRanges]{grglist}(x, ...)}.
  }
}

\value{
  An \link[IRanges]{RleList} object \emph{parallel} to \code{transcripts},
  that is, the i-th element in it is an integer-\link[S4Vectors]{Rle}
  representing the coverage of the i-th transcript in \code{transcripts}.
  Its \code{lengths()} is guaranteed to be identical to
  \code{sum(width(transcripts))}. The names and metadata columns on
  \code{transcripts} are propagated to it.
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \code{\link{extractTranscriptSeqs}} for extracting transcript
          (or CDS) sequences from chromosome sequences.

    \item \code{\link{transcriptLengths}} for extracting the transcript
          lengths from a \link{TxDb} object.

    \item The \link[IRanges]{RleList} class defined and documented in the
          \pkg{IRanges} package.

    \item The \link[GenomicRanges]{GRangesList} class defined and documented
          in the \pkg{GenomicRanges} package.

    \item The \code{\link[GenomicRanges]{coverage}} methods defined in the
          \pkg{GenomicRanges} package.

    \item The \code{\link{exonsBy}} function for extracting exon ranges
          grouped by transcript.

    \item \code{\link[GenomicAlignments]{findCompatibleOverlaps}} in the
          \pkg{GenomicAlignments} package for finding which reads are
          \emph{compatible} with the splicing of which transcript.
  }
}

\examples{
## ---------------------------------------------------------------------
## 1. A SIMPLE EXAMPLE WITH ONLY ONE TRANSCRIPT
## ---------------------------------------------------------------------

## Get some transcripts:
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
dm3_transcripts <- exonsBy(txdb, by="tx", use.names=TRUE)
dm3_transcripts

## Let's pick up the 1st transcript: FBtr0300689. It as 2 exons and 1
## intron:
my_transcript <- dm3_transcripts["FBtr0300689"]

## Let's create 3 artifical aligned reads. We represent them as a
## GRanges object of length 3 that contains the genomic positions of
## the 3 reads. Note that these reads are simple alignments i.e. each
## of them can be represented with a single range. This would not be
## the case if they were junction reads.
my_reads <- GRanges(c("chr2L:7531-7630",
                      "chr2L:8101-8200",
                      "chr2L:8141-8240"))

## The coverage of the 3 reads on the reference genome is:
coverage(my_reads)

## As you can see, all the genomic positions in the 3 ranges participate
## to the coverage. This can be confirmed by comparing:
sum(coverage(my_reads))
## with:
sum(width(my_reads))
## They should always be the same.

## When computing the coverage on a transcript, only the part of the
## read that overlaps with the transcript participates to the coverage.
## Let's look at the individual coverage of each read on transcript
## FBtr0300689:

## The 1st read is fully contained within the 1st exon:
coverageByTranscript(my_reads[1], my_transcript)

## The 2nd and 3rd reads overlap the 2 exons and the intron:
coverageByTranscript(my_reads[2], my_transcript)
coverageByTranscript(my_reads[3], my_transcript)

## The coverage of the 3 reads together is:
coverageByTranscript(my_reads, my_transcript)

## Note that this is the sum of the individual coverages. This can be
## checked with:
stopifnot(all(
  coverageByTranscript(my_reads, my_transcript)
  ==
  Reduce("+", lapply(my_reads, coverageByTranscript, my_transcript), 0L)
))

## ---------------------------------------------------------------------
## 2. COMPUTE THE FULL TRANSCRIPTOME COVERAGE OF A SET OF ALIGNED READS
## ---------------------------------------------------------------------

## Load the aligned reads:
library(pasillaBamSubset)
library(GenomicAlignments)
reads <- readGAlignments(untreated1_chr4())

## Compute the full transcriptome coverage by calling
## coverageByTranscript() on 'dm3_transcripts':
tx_cvg <- coverageByTranscript(reads, dm3_transcripts, ignore.strand=TRUE)
tx_cvg

## A sanity check:
stopifnot(identical(lengths(tx_cvg), sum(width(dm3_transcripts))))

## We can also use pcoverageByTranscript() to compute 'tx_cvg'.
## For this we first create a GAlignmentsList object "parallel" to
## 'dm3_transcripts' where the i-th list element contains the aligned
## reads that overlap with the i-th transcript:
hits <- findOverlaps(reads, dm3_transcripts, ignore.strand=TRUE)
tx2reads <- setNames(as(t(hits), "List"), names(dm3_transcripts))
reads_by_tx <- extractList(reads, tx2reads)  # GAlignmentsList object
reads_by_tx

## Call pcoverageByTranscript():
tx_cvg2 <- pcoverageByTranscript(reads_by_tx, dm3_transcripts,
                                 ignore.strand=TRUE)
stopifnot(identical(tx_cvg, tx_cvg2))

## A more meaningful coverage is obtained by counting for each
## transcript only the reads that are *compatible* with its splicing:
compat_hits <- findCompatibleOverlaps(reads, dm3_transcripts)
tx2reads <- setNames(as(t(compat_hits), "List"), names(dm3_transcripts))
compat_reads_by_tx <- extractList(reads, tx2reads)

tx_compat_cvg <- pcoverageByTranscript(compat_reads_by_tx,
                                       dm3_transcripts,
                                       ignore.strand=TRUE)
## A sanity check:
stopifnot(all(all(tx_compat_cvg <= tx_cvg)))

## ---------------------------------------------------------------------
## 3. COMPUTE CDS COVERAGE OF A SET OF ALIGNED READS
## ---------------------------------------------------------------------

## coverageByTranscript() can also be used to compute CDS coverage:
cds <- cdsBy(txdb, by="tx", use.names=TRUE)
cds_cvg <- coverageByTranscript(reads, cds, ignore.strand=TRUE)
cds_cvg

## A sanity check:
stopifnot(identical(lengths(cds_cvg), sum(width(cds))))

## ---------------------------------------------------------------------
## 4. ALTERNATIVELY, THE CDS COVERAGE CAN BE OBTAINED FROM THE
##    TRANSCRIPT COVERAGE BY TRIMMING THE 5' AND 3' UTRS
## ---------------------------------------------------------------------

tx_lens <- transcriptLengths(txdb, with.utr5_len=TRUE, with.utr3_len=TRUE)
stopifnot(identical(tx_lens$tx_name, names(tx_cvg)))  # sanity

## Keep the rows in 'tx_lens' that correspond to a list element in
## 'cds_cvg' and put them in the same order as in 'cds_cvg':
m <- match(names(cds_cvg), names(tx_cvg))
tx_lens <- tx_lens[m, ]
utr5_width <- tx_lens$utr5_len
utr3_width <- tx_lens$utr3_len

trimListElements <- function(x, ltrim=0, rtrim=0)
{
    x_lens <- lengths(x)
    n1 <- pmax(x_lens - rtrim, 0)
    n2 <- pmax(n1 - ltrim, 0)
    ptail(phead(x, n=n1), n=n2)
}

cds_cvg2 <- trimListElements(tx_cvg[m], utr5_width, utr3_width)

## A sanity check:
stopifnot(identical(cds_cvg2, cds_cvg))
}

\keyword{manip}
