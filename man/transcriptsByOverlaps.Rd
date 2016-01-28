\name{transcriptsByOverlaps}

\alias{transcriptsByOverlaps}
\alias{transcriptsByOverlaps,TxDb-method}
\alias{exonsByOverlaps}
\alias{exonsByOverlaps,TxDb-method}
\alias{cdsByOverlaps}
\alias{cdsByOverlaps,TxDb-method}

\title{
  Extract genomic features from an object based on their by genomic location
}
\description{
  Generic functions to extract genomic features for specified genomic
  locations.
  This page documents the methods for \link{TxDb} objects only.
}
\usage{
transcriptsByOverlaps(x, ranges,
                      maxgap = 0L, minoverlap = 1L,
                      type = c("any", "start", "end"), ...)
\S4method{transcriptsByOverlaps}{TxDb}(x, ranges,
                      maxgap = 0L, minoverlap = 1L,
                      type = c("any", "start", "end"),
                      columns = c("tx_id", "tx_name"))

exonsByOverlaps(x, ranges,
                maxgap = 0L, minoverlap = 1L,
                type = c("any", "start", "end"), ...)
\S4method{exonsByOverlaps}{TxDb}(x, ranges,
                maxgap = 0L, minoverlap = 1L,
                type = c("any", "start", "end"),
                columns = "exon_id")

cdsByOverlaps(x, ranges,
              maxgap = 0L, minoverlap = 1L,
              type = c("any", "start", "end"), ...)
\S4method{cdsByOverlaps}{TxDb}(x, ranges,
              maxgap = 0L, minoverlap = 1L,
              type = c("any", "start", "end"),
              columns = "cds_id")
}
\arguments{  
  \item{x}{A \link{TxDb} object.}
  \item{...}{Arguments to be passed to or from methods.}
  \item{ranges}{A \link[GenomicRanges]{GRanges} object to restrict the output.}
  \item{type}{How to perform the interval overlap operations of the
    \code{ranges}. See the
    \code{\link[GenomicRanges:findOverlaps-methods]{findOverlaps}} manual page
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
  using the results of \code{\link[IRanges]{findOverlaps}}
  calls based on the specified \code{ranges}.
}
\value{ a GRanges object }
\author{
  P. Aboyoun
}
\seealso{
  \itemize{
    \item \code{\link{transcripts}} and \code{\link{transcriptsBy}}
          for more ways to extract genomic features
          from a \link{TxDb} object.
    \item \link[GenomicFeatures]{select-methods} for how to use the
          simple "select" interface to extract information from a
          \link{TxDb} object.
    \item \code{\link{id2name}} for mapping \link{TxDb} internal ids
          to external names for a given feature type.
    \item The \link{TxDb} class.
  }
}
\examples{
  txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite",
                                   package="GenomicFeatures"))
  gr <- GRanges(seqnames = rep("chr1",2),
                ranges = IRanges(start=c(500,10500), end=c(10000,30000)),
                strand = strand(rep("-",2)))
  transcriptsByOverlaps(txdb, gr)
}

\keyword{methods}
