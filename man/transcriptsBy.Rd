\name{transcriptsBy}

\alias{transcriptsBy}
\alias{exonsBy}
\alias{cdsBy}
\alias{fiveUTRsByTranscript}
\alias{threeUTRsByTranscript}

\title{
  Extract and group genomic features of a given type
}
\description{
  These functions extract genomic features of a given type
  from a \link{TranscriptDb} object and group them according to
  a user-specified criteria.
}
\usage{
  transcriptsBy(txdb, by = c("gene", "exon", "cds"))
  exonsBy(txdb, by = c("tx", "gene"))
  cdsBy(txdb, by = c("tx", "gene"))
  fiveUTRsByTranscript(txdb)
  threeUTRsByTranscript(txdb)
}
\arguments{
  \item{txdb}{A \link{TranscriptDb} object.}
  \item{by}{One of \code{"gene"}, \code{"exon"}, \code{"cds"} or \code{"tx"}.
    Determines the grouping.}
}
\details{
  These functions return a \link[GenomicRanges]{GRangesList} object
  where the ranges within each of the elements are ordered according
  to the following rule:

  When using \code{exonsBy} and \code{cdsBy} with \code{by = "tx"},
  the ranges are returned in the order they appear in the transcript,
  i.e. order by the splicing.exon_rank field in \code{txdb}'s internal
  database. In all other cases, the ranges will be ordered by chromosome,
  strand, start, and end values.
}
\value{A \link[GenomicRanges]{GRangesList} object.}
\author{
  M. Carlson, P. Aboyoun and H. Pages
}
\seealso{
  \link{TranscriptDb}, \code{\link{transcripts}},
  \code{\link{transcriptsByOverlaps}}
}
\examples{
  txdb_file <- system.file("extdata", "UCSC_knownGene_sample.sqlite",
                           package="GenomicFeatures")
  txdb <- loadFeatures(txdb_file)

  ## Get the transcript ranges grouped by gene
  transcriptsBy(txdb, "gene")

  ## Get the exons grouped by gene
  exonsBy(txdb, "gene")

  ## Get the cds ranges grouped by transcript
  cdsBy(txdb, "tx")

  ## Get the 5' UTRs grouped by transcript
  fiveUTRsByTranscript(txdb)
}
