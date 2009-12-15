
\name{getTranscripts}
\alias{getTranscripts}
\alias{getExons}
\title{
  Functions to retrieve genomic features for ranges of interest.
}
\description{
  Functions that retrieve genomic features of interest from a
  transcripts database.
}
\usage{
getTranscripts(transcript,ranges,chromosome,strand)
getExons(transcript,ranges,chromosome,strand)
}
\arguments{
  \item{transcript}{a transcripts object, which is a connection to a DB.}
  \item{ranges}{An IRanges object to specify which ranges to search for
annotations}
  \item{chromosome}{the chromosome to search on 'chr1'.}
  \item{strand}{the strand to search on +/1}
}
\details{These functions return relevant annotation results in a data
frame based on the regions searched.}
\value{ All of these functions return a data frame. }
\author{
  M. Carlson
}

\examples{
\dontrun{
  tx <- UCSCTranscripts(organism="human")
  foo <- IRanges(start=c(500,10500), end=c(10000,30000))
  getTranscripts(tx,foo,"chr1","-")
}
}
