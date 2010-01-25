\name{transcripts}

\alias{transcripts}
\alias{exons}
\alias{cds}

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
  cds(txdb, vals)
}
\arguments{  
  \item{txdb}{A \link{TranscriptDb} object.}
  \item{vals}{Either NULL or a named list of vectors to be used to restrict
    the output.
    For \code{transcripts} the element names can be "gene_id", "tx_id",
    "tx_name", "tx_chrom", and "tx_strand".
    For \code{exons} the element names can be "exon_id", "exon_chrom", and
    "exon_strand".
    For \code{cds} the element names can be "cds_id", "cds_chrom", and
    "cds_strand".}
  \item{columns}{columns to include in the output. Can be NULL or
    a combination of "tx_id", "tx_name", "gene_id", "exon_id", and
    "cds_id".}
}
\details{
  These are the main functions for retrieving transcript information from
  a \link{TranscriptDb} object. They can restrict the output based on
  categorical information. To restrict the output based on interval
  information, use the \code{\link{transcriptsByRanges}},
  \code{\link{exonsByRanges}}, and \code{\link{cdsByRanges}}
  functions.
}
\value{ a RangedData object }
\author{
  M. Carlson and P. Aboyoun
}
\seealso{
  \link{TranscriptDb}, \code{\link{transcriptsByRanges}},
  \code{\link{bindTranscripts}}
}
\examples{
  txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                   package="GenomicFeatures"))
  vals <- list(tx_chrom = c("chr3", "chr5"), tx_strand = "+")
  transcripts(txdb, vals)
}
