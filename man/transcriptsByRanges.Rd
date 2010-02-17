\name{transcriptsByRanges}

\alias{transcriptsByRanges}
\alias{exonsByRanges}
\alias{cdsByRanges}

\title{
  Retrieving genomic features from a TranscriptDb object using a restriction
  by genomic location.
}
\description{
  Functions to retrieve genomic features for specified genomic locations.
}
\usage{
  transcriptsByRanges(txdb, ranges,
                      restrict = c("any", "start", "end", "within", "equal"),
                      columns = c("tx_id","tx_name"))
  exonsByRanges(txdb, ranges,
                restrict = c("any", "start", "end", "within", "equal"))
  cdsByRanges(txdb, ranges,
              restrict = c("any", "start", "end", "within", "equal"))
}
\arguments{  
  \item{txdb}{A \link{TranscriptDb} object.}
  \item{ranges}{A \link[BSgenome]{GenomicFeature} object to restrict the output.}
  \item{restrict}{How to perform the interval overlap operations of the
    \code{ranges}. See the
    \code{\link[IRanges:IntervalTree-class]{findOverlaps}} manual page
    in the IRanges package for more information.}
  \item{columns}{columns to include in the output. Can be NULL or
    a combination of "tx_id", "tx_name", "gene_id", "exon_id", and
    "cds_id".}
}
\details{
  These functions subset the results of \code{\link{transcripts}},
  \code{\link{exons}}, and \code{\link{cds}} function calls with
  using the results of \code{\link[IRanges:IntervalTree-class]{findOverlaps}}
  calls based on the specified \code{ranges}.
}
\value{ a GenomicFeature object }
\author{
  P. Aboyoun and M. Carlson
}
\seealso{
  \link{TranscriptDb}, \code{\link{transcripts}}, \code{\link{bindTranscripts}}
}
\examples{
  library(IRanges)
  library(BSgenome)
  txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                   package="GenomicFeatures"))
  gf <- GenomicFeature(seqnames = rep("chr1",2),
                       ranges = IRanges(start=c(500,10500), end=c(10000,30000)),
                       strand = strand(rep("-",2)))
  transcriptsByRanges(txdb, gf)
}
