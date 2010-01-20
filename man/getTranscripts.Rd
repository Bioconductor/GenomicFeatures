\name{transcriptsByRanges}

\alias{transcriptsByRanges}
\alias{transcriptsByRanges,missing-method}
\alias{transcriptsByRanges,RangedData-method}
\alias{exonsByRanges}
\alias{exonsByRanges,missing-method}
\alias{exonsByRanges,RangedData-method}

\title{
  Retrieving genomic features from a TranscriptDb object.
}
\description{
  Methods to retrieve genomic features for ranges of interest from a
  TranscriptDb object.
}
\usage{
  transcriptsByRanges(txdb, ranges=NULL, restrict, columns)
  exonsByRanges(txdb, ranges=NULL, restrict, columns)
}
\arguments{  
  \item{txdb}{A \link{TranscriptDb} object.}
  \item{ranges}{A RangedData object to specify which ranges to search
    for annotations.}
  \item{restrict}{How to test the ranges.  Can be "any", "start",
    "during", "finish" or "equal" to indicate the kind of overlap
    allowed between the ranges to be mapped and the ranges in the DB.
    The default value is "any".}
  \item{columns}{Which data columns to return with the rangedData object.}
}
\details{These functions return relevant annotation results in a Ranged
Data object based on the regions searched.  The \code{exonsByRanges} methods are
for extracting data based on which exons are contained within a set of
ranges, while a call to \code{transcriptsByRanges} will extract data based
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
  transcriptsByRanges(txdb,rd)
}
}
