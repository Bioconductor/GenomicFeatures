\name{extractTranscriptCoverage}

\alias{extractTranscriptCoverage}

\title{Extract transcript or CDS coverage from a set of aligned reads}

\description{
  \code{extractTranscriptCoverage} extracts transcript or CDS coverage from
  a set of aligned reads.
}

\usage{
extractTranscriptCoverage(reads, transcripts)
}

\arguments{
  \item{reads}{
    An object representing a set of aligned reads. Can be a
    \link[Rsamtools]{BamFile}, \link[GenomicAlignments]{GAlignments},
    \link[GenomicAlignments]{GAlignmentPairs},
    \link[GenomicAlignments]{GAlignmentsList}, or
    \link[GenomicRanges]{GRangesList} object, or, more generally, any
    object for which \code{\link[GenomeInfoDb]{seqinfo}()} and
    \code{\link[GenomicRanges]{coverage}()} are supported (note that,
    for such objects, \code{coverage()} is expected to return an
    \link[IRanges]{RleList} object whose names are \code{seqlevels(reads)}).
  }
  \item{transcripts}{
    A \link[GenomicRanges]{GRangesList} object representing the exons of
    each transcript for which to extract coverage. For each transcript, the
    exons must be ordered by \emph{ascending rank}, that is, by their position
    in the transcript. This means that, for a transcript located on the minus
    strand, the exons should typically be ordered by descending position on
    the reference genome. If \code{transcripts} was obtained with
    \code{\link{exonsBy}}, then the exons are guaranteed to be ordered by
    ascending rank. See \code{?\link{exonsBy}} for more information.

    Alternatively \code{transcripts} can be any object for which
    \code{\link{exonsBy}} is implemented (e.g. a \link{TxDb} object). In
    that case, it's replaced by the \link[GenomicRanges]{GRangesList} object
    returned by \code{\link{exonsBy}(transcripts, by="tx", use.names=TRUE)}.
  }
}

\value{
  An \link[IRanges]{RleList} object \emph{parallel} to \code{transcripts},
  that is, the i-th element in it is an integer-\link[S4Vectors]{Rle}
  representing the coverage of the i-th transcript in \code{transcripts}.
  Its \code{elementLengths()} is guaranteed to be identical to
  \code{sum(width(transcripts))}. The names and metadata columns on
  \code{transcripts} are propagated to it.
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item \code{\link{extractTranscriptSeqs}} for extracting transcript or
          CDS sequences from chromosome sequences.

    \item The \link[IRanges]{RleList} class defined and documented in the
          \pkg{IRanges} package.

    \item The \link[GenomicRanges]{GRangesList} class defined and documented
          in the \pkg{GenomicRanges} package.

    \item The \code{\link[GenomicRanges]{coverage}} methods defined in the
          \pkg{GenomicRanges} package.

    \item The \code{\link{exonsBy}} function for extracting exon ranges
          grouped by transcript.

    \item The \link{TxDb} class.
  }
}

\examples{
## Load the aligned reads:
library(pasillaBamSubset)
library(GenomicAlignments)
reads <- readGAlignments(untreated1_chr4())

## Load the transcripts:
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
transcripts <- exonsBy(txdb, by="tx", use.names=TRUE)

## Extract the transcript coverage:
cvg_by_tx <- extractTranscriptCoverage(reads, transcripts)
cvg_by_tx

## A sanity check:
stopifnot(identical(elementLengths(cvg_by_tx), sum(width(transcripts))))
}

\keyword{manip}
