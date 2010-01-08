\name{getTranscripts}

\alias{getTranscripts}
\alias{getTranscripts,TranscriptDb,missing-method}
\alias{getTranscripts,TranscriptDb,IRanges-method}
\alias{getTranscripts,TranscriptDb,RangesList-method}
\alias{getTranscripts,TranscriptDb,RangedData-method}
\alias{getExons}
\alias{getExons,TranscriptDb,missing-method}
\alias{getExons,TranscriptDb,IRanges-method}
\alias{getExons,TranscriptDb,RangesList-method}
\alias{getExons,TranscriptDb,RangedData-method}

\title{
  Retrieveing genomic features from a TranscriptDb object.
}
\description{
  Methods to retrieve genomic features for ranges of interest from a
  TranscriptDb object.
}
\usage{
  getTranscripts(ann, ranges=NULL, chromosome=NULL, strand=NULL,
                 rangeRestr="either", expand=FALSE)
  getExons(ann, ranges=NULL, chromosome=NULL, strand=NULL,
           rangeRestr="either", expand=FALSE)
}
\arguments{
  \item{ann}{A \link{TranscriptDb} object.}
  \item{ranges}{An IRanges, RangesList or RangedData object to specify
    which ranges to search for annotations.}
  \item{chromosome}{the chromosome to search on 'chr1'.}
  \item{strand}{the strand to search on +/1}
  \item{rangeRestr}{How to test the ranges.  Either "both", "start",
    "end" or "either" of any range can be required to land within an
    annotation.}
  \item{expand}{do additional joins to include transcript information or
    exon information depending on whether \code{getTranscripts} or
    \code{getExons} is called.}
}
\details{These functions return relevant annotation results in a data
frame based on the regions searched.  The \code{getExons} methods are
for extracting data based on which exons are contained within a set of
ranges, while a call to \code{getTranscripts} will extract data based
on the transcript boundaries being contained in the ranges.
The chromosome and strand information parameters may be optional if
that information is already provided by the ranges object supplied.}
\value{ All of these functions return a data frame. }
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
  getTranscripts(txdb, rd)
}
}
