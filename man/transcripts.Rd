\name{transcripts}

\alias{transcripts}
\alias{exons}
\alias{cds}

\title{
  Retrieving genomic features from a TranscriptDb object.
}
\description{
  Functions to retrieve genomic features from a \link{TranscriptDb} object.
}
\usage{
transcripts(txdb, vals=NULL, columns=c("tx_id", "tx_name"))
exons(txdb, vals=NULL, columns="exon_id")
cds(txdb, vals=NULL)
}
\arguments{  
  \item{txdb}{
    A \link{TranscriptDb} object.
  }
  \item{vals}{
    Either \code{NULL} or a named list of vectors to be used to
    restrict the output.
    For \code{transcripts} the element names can be \code{"gene_id"},
    \code{"tx_id"}, \code{"tx_name"}, \code{"tx_chrom"},
    and \code{"tx_strand"}.
    For \code{exons} the element names can be \code{"exon_id"},
    \code{"exon_chrom"}, and \code{"exon_strand"}.
    For \code{cds} the element names can be \code{"cds_id"},
    \code{"cds_chrom"}, and \code{"cds_strand"}.
  }
  \item{columns}{
    Columns to include in the output.
    Must be \code{NULL} or a character vector with values in:
    \itemize{
    \item \code{"tx_id"}, \code{"tx_name"}, \code{"gene_id"},
          \code{"exon_id"}, and \code{"cds_id"} for \code{transcripts}.
    \item \code{"gene_id"},
          \code{"tx_id"}, \code{"tx_name"},
          \code{"tx_chrom"}, \code{"tx_strand"},
          \code{"exon_id"}, \code{"exon_name"},
          \code{"cds_id"}, \code{"cds_name"},
          \code{"cds_chrom"} and \code{"cds_strand"} for \code{exons}.
    }
  }
}
\details{
  These are the main functions for retrieving transcript information from
  a \link{TranscriptDb} object. They can restrict the output based on
  categorical information. To restrict the output based on interval
  information, use the \code{\link{transcriptsByOverlaps}},
  \code{\link{exonsByOverlaps}}, and \code{\link{cdsByOverlaps}}
  functions.
}
\value{ a GRanges object }
\author{
  M. Carlson, P. Aboyoun and H. Pages
}
\seealso{
  \link{TranscriptDb}, \code{\link{transcriptsByOverlaps}}
}
\examples{
  txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                   package="GenomicFeatures"))
  vals <- list(tx_chrom = c("chr3", "chr5"), tx_strand = "+")
  transcripts(txdb, vals)
}
