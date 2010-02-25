\name{transcriptsBy}

\alias{transcriptsBy}
\alias{exonsBy}
\alias{cdsBy}

\title{
  Retrieve a GRangeList from a TranscriptDb object grouped by gene,
  transcript, exon or cds.
}
\description{
  Functions to retrieve a GRangesList from a \link{TranscriptDb}
  object grouped by gene, transcript, exon or cds.
}
\usage{
  transcriptsBy(txdb, groupBy)
  exonsBy(txdb, groupBy)
  cdsBy(txdb, groupBy)
}
\arguments{  
  \item{txdb}{A \link{TranscriptDb} object.}
  \item{groupBy}{Any of "tx", "gene", "exon" or "cds".  Used to
  determine how the results that are returned are grouped.
}
}
\details{
  These are functions for retrieving grouped transcript information from
  a \link{TranscriptDb} object.  A GRangesList class is designed to hold
  collections of GRanges objects and thus group things like transcripts
  (which typically have several range elements each) in a sensible way.
  These accessors return annotations grouped in the way that the user
  requests and packaged into a convenient GRangesList object
}
\value{ a GRangesList object }
\author{
  M. Carlson
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
