\name{getTranscripts}

\alias{getTranscripts}
\alias{getTranscripts,missing-method}
\alias{getTranscripts,RangedData-method}
\alias{getExons}
\alias{getExons,missing-method}
\alias{getExons,RangedData-method}

\title{
  Retrieving genomic features from a TranscriptDb object.
}
\description{
  Methods to retrieve genomic features for ranges of interest from a
  TranscriptDb object.
}
\usage{
  getTranscripts(ranges=NULL, ann, rangeRestr="either", expand=FALSE)
  getExons(ranges=NULL, ann, rangeRestr="either", expand=FALSE)
}
\arguments{
  \item{ranges}{A RangedData object to specify which ranges to search for annotations.}
  \item{ann}{A \link{TranscriptDb} object.}
  \item{rangeRestr}{How to test the ranges.  Either "both", "start",
    "end" or "either" of any range can be required to land within an
    annotation.}
  \item{expand}{do additional joins to include transcript information or
    exon information depending on whether \code{getTranscripts} or
    \code{getExons} is called.}
}
\details{These functions return relevant annotation results in a Ranged
Data object based on the regions searched.  The \code{getExons} methods are
for extracting data based on which exons are contained within a set of
ranges, while a call to \code{getTranscripts} will extract data based
on the transcript boundaries being contained in the ranges.}
\value{ All of these methods should return a RangedData object. }
\author{
  M. Carlson
}
\seealso{
  \link{TranscriptDb}
}
\examples{
\dontrun{
  txdb <- loadFeatures(system.file("extdata", "HG18test.sqlite",
                                   package="GenomicFeatures"))
  library(IRanges)
  rd <- RangedData(ranges = IRanges(start=c(500,10500),end=c(10000,30000) ),
                   space = rep("chr1",2), 
                   strand = rep("-",2))
  getTranscripts(rd, txdb)
}
}
