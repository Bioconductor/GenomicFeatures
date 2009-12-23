\name{getTranscripts}

\alias{getTranscripts}
\alias{getTranscripts,TranscriptAnnotation,IRanges-method}
\alias{getTranscripts,TranscriptAnnotation,RangesList-method}
\alias{getTranscripts,TranscriptAnnotation,RangedData-method}
\alias{getExons}
\alias{getExons,TranscriptAnnotation,IRanges-method}
\alias{getExons,TranscriptAnnotation,RangesList-method}
\alias{getExons,TranscriptAnnotation,RangedData-method}

\title{
  Retrieveing genomic features from a TranscriptAnnotation object.
}
\description{
  Methods to retrieve genomic features for ranges of interest from a
  TranscriptAnnotation object.
}
\usage{
  getTranscripts(ann, ranges=NULL, chromosome=NULL, strand=NULL,
                 rangeRestr="either", expand=FALSE, showSQL=FALSE)
  getExons(ann, ranges=NULL, chromosome=NULL, strand=NULL,
           rangeRestr="either", expand=FALSE, showSQL=FALSE)
}
\arguments{
  \item{ann}{A \link{TranscriptAnnotation} object.}
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
  \item{showSQL}{Show the SQL query sent to the DB?}
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
  \link{TranscriptAnnotation}
}
\examples{
\dontrun{
  txann <- UCSCTranscripts()
  foo <- IRanges(start=c(500,10500), end=c(10000,30000))
  getTranscripts(txann, foo, "chr1","-")
}
}
