\name{transcriptsBy}

\alias{transcriptsBy}
\alias{exonsBy}
\alias{cdsBy}

\title{
  Retrieving genomic features from a TranscriptDb object that are
  grouped by gene, transcript, exon or cds IDs.
}
\description{
  Functions to retrieve genomic features from a \link{TranscriptDb}
  that are grouped by gene, transcript, exon or cds.
}
\usage{
  transcriptsBy(txdb, groupBy = c("gene", "exon", "cds"))
  exonsBy(txdb, groupBy = c("tx", "gene"))
  cdsBy(txdb, groupBy = c("tx", "gene"))
}
\arguments{
  \item{txdb}{A \link{TranscriptDb} object.}
  \item{groupBy}{One of "tx", "gene", "exon" or "cds". Determines
    how to group the results.}
}
\details{
  These function return a GRangesList object where the ranges within
  each of the elements are ordered according to the following rule:

  When using \code{exonsBy} and \code{cdsBy} with \code{groupBy = "tx"},
  the ranges are returned in the order they appear in the transcript.
  In all other cases, the ranges will be ordered by chromosome, start,
  and end values.
}
\value{ a GRangesList object }
\author{
  M. Carlson and P. Aboyoun
}
\seealso{
  \link{TranscriptDb}, \code{\link{transcripts}},
  \code{\link{transcriptsByRanges}}, \code{\link{bindTranscripts}}
}
\examples{
  txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                   package="GenomicFeatures"))

  ## Get the transcript ranges grouped by gene
  transcriptsBy(txdb, "gene")

  ## Get the exons grouped by gene
  exonsBy(txdb, "gene")

  ## Get the cds ranges grouped by transcript
  cdsBy(txdb, "tx")
}
