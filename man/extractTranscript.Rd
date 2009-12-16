\name{getTranscripts}
\alias{getTranscripts}
\alias{getTranscripts,IRanges-method}
\alias{getTranscripts,RangedData-method}
\alias{getTranscripts,RangesList-method}
\alias{getExons}
\alias{getExons,IRanges-method}
\alias{getExons,RangedData-method}
\alias{getExons,RangesList-method}
\title{
  Methods to retrieve genomic features for ranges of interest.
}
\description{
  Methods that retrieve genomic features of interest from a
  transcripts database.
}
\usage{
  getTranscripts(x, transcript, chromosome, strand, rangeRestr, expand,
                 showSQL)
  getExons(x, transcript, chromosome, strand, rangeRestr, expand,
           showSQL)
}
\arguments{
  \item{x}{An IRanges object to specify which ranges to search for
    annotations}
  \item{transcript}{a transcripts object, which is a connection to a DB.}
  \item{chromosome}{the chromosome to search on 'chr1'.}
  \item{strand}{the strand to search on +/1}
  \item{rangeRestr}{How to test the ranges.  Either "both", "start",
    "end" or "either" of any range can be required to land within an
    annotation.}
  \item{expand}{do additional joins to include transcript information or
    exon information depending on whether getTranscripts or getExons is
    called.}
  \item{showSQL}{Show the SQL sent to the DB?}
}
\details{These functions return relevant annotation results in a data
frame based on the regions searched.  The getExons() methods are for
extracting data based on which exons are contained within a set of
ranges, while a call to getTranscripts() will extract data based on the
transcript boundaries being contained in the ranges.  The Chromosome and
strand information parameters may be optional if that information is
already provided by the ranges object supplied.}
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
