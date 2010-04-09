\name{transcriptsByOverlaps}

\alias{transcriptsByOverlaps}
\alias{exonsByOverlaps}
\alias{cdsByOverlaps}

\title{
  Retrieving genomic features from a TranscriptDb object using a restriction
  by genomic location.
}
\description{
  Functions to retrieve genomic features for specified genomic locations.
}
\usage{
  transcriptsByOverlaps(txdb, ranges, maxgap = 0L, minoverlap = 1L,
                        type = c("any", "start", "end"),
                        columns = c("tx_id", "tx_name"))
  exonsByOverlaps(txdb, ranges, maxgap = 0L, minoverlap = 1L,
                  type = c("any", "start", "end"))
  cdsByOverlaps(txdb, ranges, maxgap = 0L, minoverlap = 1L,
                type = c("any", "start", "end"))
}
\arguments{  
  \item{txdb}{A \link{TranscriptDb} object.}
  \item{ranges}{A \link[GenomicRanges]{GRanges} object to restrict the output.}
  \item{type}{How to perform the interval overlap operations of the
    \code{ranges}. See the
    \code{\link[GenomicRanges:GRanges-findOverlaps]{findOverlaps}} manual page
    in the GRanges package for more information.}
  \item{maxgap}{A non-negative integer representing the maximum distance
    between a query interval and a subject interval.}
  \item{minoverlap}{Ignored.}
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
\value{ a GRanges object }
\author{
  P. Aboyoun
}
\seealso{
  \link{TranscriptDb}, \code{\link{transcripts}}, \code{\link{bindTranscripts}}
}
\examples{
  txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                   package="GenomicFeatures"))
  gr <- GRanges(seqnames = rep("chr1",2),
                ranges = IRanges(start=c(500,10500), end=c(10000,30000)),
                strand = strand(rep("-",2)))
  transcriptsByOverlaps(txdb, gr)
}
