\name{transcripts}

\alias{transcripts}
\alias{exons}

\title{
  Retrieving genomic features from a TranscriptDb object.
}
\description{
  Functions to retrieve genomic features for IDs of types specified by the
  vals parameter.
}
\usage{
  transcripts(txdb, vals, columns=c("tx_id", "tx_name"))
  exons(txdb, vals)
}
\arguments{  
  \item{txdb}{A \link{TranscriptDb} object.}
  \item{vals}{A named list, where the names indicate the type
    expected.  The names can be: "gene_id", "tx_id", "tx_name",
    "tx_chrom", or "tx_strand".}
  \item{columns}{Which data columns to return with the rangedData object.}
}
\details{These functions return relevant annotation results in a Ranged
Data object based on the regions searched.  The \code{exonsByRanges}
functions are for extracting data based on which exons are contained
within a set of ranges, while a call to \code{transcriptsByRanges} will
extract data based on the transcript boundaries being contained in the
ranges.}
\value{ All of these functions should return a RangedData object. }
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
  vals <- list(tx_chrom = c("chr3", "chr5"), tx_strand = "+")
  transcripts(txdb, vals)
}
}
