\name{transcriptsByOverlaps}

\alias{transcriptsByOverlaps}
\alias{transcriptsByOverlaps,TranscriptDb-method}
\alias{exonsByOverlaps}
\alias{exonsByOverlaps,TranscriptDb-method}
\alias{cdsByOverlaps}
\alias{cdsByOverlaps,TranscriptDb-method}

\title{
  Retrieving genomic features from an object using a restriction
  by genomic location.
}
\description{
  Generic functions to retrieve genomic features for specified genomic
  locations.
  This page documents the methods for \link{TranscriptDb} objects only.
}
\usage{
transcriptsByOverlaps(x, ranges,
                      maxgap = 0L, minoverlap = 1L,
                      type = c("any", "start", "end"), ...)
\S4method{transcriptsByOverlaps}{TranscriptDb}(x, ranges,
                      maxgap = 0L, minoverlap = 1L,
                      type = c("any", "start", "end"),
                      columns = c("tx_id", "tx_name"))

exonsByOverlaps(x, ranges,
                maxgap = 0L, minoverlap = 1L,
                type = c("any", "start", "end"), ...)
\S4method{exonsByOverlaps}{TranscriptDb}(x, ranges,
                maxgap = 0L, minoverlap = 1L,
                type = c("any", "start", "end"),
                columns = "exon_id")

cdsByOverlaps(x, ranges,
              maxgap = 0L, minoverlap = 1L,
              type = c("any", "start", "end"), ...)
\S4method{cdsByOverlaps}{TranscriptDb}(x, ranges,
              maxgap = 0L, minoverlap = 1L,
              type = c("any", "start", "end"),
              columns = "cds_id")
}
\arguments{  
  \item{x}{A \link{TranscriptDb} object.}
  \item{...}{Arguments to be passed to or from methods.}
  \item{ranges}{A \link[GenomicRanges]{GRanges} object to restrict the output.}
  \item{type}{How to perform the interval overlap operations of the
    \code{ranges}. See the
    \code{\link[GenomicRanges:GRanges-findOverlaps]{findOverlaps}} manual page
    in the GRanges package for more information.}
  \item{maxgap}{A non-negative integer representing the maximum distance
    between a query interval and a subject interval.}
  \item{minoverlap}{Ignored.}
  \item{columns}{Columns to include in the output.
    See \code{?\link{transcripts}} for the possible values.}
}
\details{
  These functions subset the results of \code{\link{transcripts}},
  \code{\link{exons}}, and \code{\link{cds}} function calls with
  using the results of \code{\link[IRanges:IntervalTree-class]{findOverlaps}}
  calls based on the specified \code{ranges}.
}
\value{ a GRanges object }
\author{
  P. Aboyoun
}
\seealso{
  \link{TranscriptDb}, \code{\link{transcripts}}
}
\examples{
  txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                   package="GenomicFeatures"))
  gr <- GRanges(seqnames = rep("chr1",2),
                ranges = IRanges(start=c(500,10500), end=c(10000,30000)),
                strand = strand(rep("-",2)))
  transcriptsByOverlaps(txdb, gr)
}
