\name{bindTranscripts}

\alias{bindTranscripts}
\alias{bindTranscripts,RangedData-method}
\alias{bindExons}
\alias{bindExons,RangedData-method}

\title{
  Mapping to genomic features in a TranscriptDb object.
}
\description{
  Methods to map a set of ranges to the matching features in a
  TranscriptDb object.
}
\usage{
  bindTranscripts(txdb, ranges, restrict="any", columns=c("tx_id","tx_name"))
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
\details{These Methods return relevant annotation results for the ranges
  specified in the RangedData object that is passed in.  If the user
  wishes to restrict this search to specific chromosomes and strands,
  then that information should be included in the passed in RangedData
  object by filling in the space section for that object and also adding
  a column for the strand.  If this information is missing from the
  input RangedData object, then all possible chromosomes/strands will be
  searched.
  The object returned will also be a RangedData object that
  contains the original ranges, along with a unique transcript ID for
  getting the rest of the transcript data from the DB.  The returned
  object will always indicate the chromosome and strand that is
  associated with the newly attached transcript IDs.}
\author{
  M. Carlson
}
\seealso{
  \link{TranscriptDb}
}
\examples{
\dontrun{
  txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                   package="GenomicFeatures"))
  library(IRanges)
  rd <- RangedData(ranges = IRanges(start=c(500,10500),end=c(10000,30000)),
                   space = rep("chr1",2), 
                   strand = rep("-",2))
  bindTranscripts(txdb, rd)
}
}
