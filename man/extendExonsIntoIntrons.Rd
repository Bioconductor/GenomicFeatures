\name{extendExonsIntoIntrons}

\alias{extendExonsIntoIntrons}

\title{
  Extend exons by a given number of bases into their adjacent introns
}

\description{
  \code{extendExonsIntoIntrons} extends the supplied exons by a given
  number of bases into their adjacent introns.
}

\usage{
extendExonsIntoIntrons(ex_by_tx, extent=2)
}

\arguments{
  \item{ex_by_tx}{
    A \link[GenomicRanges]{GRangesList} object containing exons grouped
    by transcript. This must be an object as returned by
    \code{\link{exonsBy}(txdb, by="tx")}, that is:
    \itemize{
      \item each list element in \code{ex_by_tx} must be a
            \link[GenomicRanges]{GRanges} object representing the
            exons of a given transcript;
      \item the exons in each list element must be ordered by ascending
            rank with respect to their transcript.
    }
  }
  \item{extent}{
    Size of the extent in number of bases. 2 by default.

    The first exon in a transcript will be extended by that amount on
    its 3' side only. The last exon in a transcript will be extended by
    that amount on its 5' side only. All other exons (i.e. intermediate
    exons) will be extended by that amount on \emph{each} side.

    Note that exons that belong to a single-exon transcript don't get
    extended.

    The default value of 2 corresponds to inclusion of the donor/acceptor
    intronic regions (typically GT/AG).
  }
}

\value{
  A copy of \link[GenomicRanges]{GRangesList} object \code{ex_by_tx}
  where the original exon ranges have been extended.

  Names and metadata columns on \code{ex_by_tx} are propagated to the
  result.
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \code{\link{transcripts}}, \code{\link{transcriptsBy}},
          and \code{\link{transcriptsByOverlaps}}, for extracting
          genomic feature locations from a \link{TxDb}-like object.

    \item \code{\link{exonicParts}} and \code{\link{intronicParts}} for
          extracting non-overlapping exonic or intronic parts from a
          TxDb-like object.

    \item \code{\link{extractTranscriptSeqs}} for extracting transcript
          (or CDS) sequences from chromosome sequences.

    \item The \link{TxDb} class.
  }
}

\examples{
## With toy transcripts:
ex_by_tx <- GRangesList(
    TX1="chr1:10-20:+",
    TX2=c("chr1:10-20:+", "chr1:50-75:+"),
    TX3=c("chr1:10-20:+", "chr1:50-75:+", "chr1:100-120:+"),
    TX4="chr1:10-20:-",
    TX5=c("chr1:10-20:-", "chr1:50-75:-"),
    TX6=c("chr1:10-20:-", "chr1:50-75:-", "chr1:100-120:-")
)

extended <- extendExonsIntoIntrons(ex_by_tx, extent=2)
extended[1:3]
extended[4:6]

## With real-world transcripts:
library(TxDb.Celegans.UCSC.ce11.ensGene)
txdb <- TxDb.Celegans.UCSC.ce11.ensGene
ex_by_tx <- exonsBy(txdb, by="tx")
ex_by_tx

extendExonsIntoIntrons(ex_by_tx, extent=2)

## Sanity check:
stopifnot(identical(extendExonsIntoIntrons(ex_by_tx, extent=0), ex_by_tx))
}

\keyword{manip}
