\name{transcripts}

\alias{transcripts}
\alias{transcripts,data.frame-method}
\alias{transcripts,TranscriptDb-method}
\alias{exons}
\alias{exons,data.frame-method}
\alias{exons,TranscriptDb-method}
\alias{cds}
\alias{cds,TranscriptDb-method}

\title{
  Retrieving genomic features from a TranscriptDb object.
}

\description{
  Generic functions to retrieve genomic features from an object.
  This page documents the methods for \link{TranscriptDb} objects only.
}

\usage{
transcripts(x, ...)
\S4method{transcripts}{TranscriptDb}(x, vals=NULL, columns=c("tx_id", "tx_name"))

exons(x, ...)
\S4method{exons}{TranscriptDb}(x, vals=NULL, columns="exon_id")

cds(x, ...)
\S4method{cds}{TranscriptDb}(x, vals=NULL, columns="cds_id")
}

\arguments{  
  \item{x}{
    A \link{TranscriptDb} object.
  }
  \item{...}{
    Arguments to be passed to or from methods.
  }
  \item{vals}{
    Either \code{NULL} or a named list of vectors to be used to
    restrict the output. Valid names for this list are:
    \code{"gene_id"},
    \code{"tx_id"}, \code{"tx_name"},
    \code{"tx_chrom"}, \code{"tx_strand"},
    \code{"exon_id"}, \code{"exon_name"},
    \code{"exon_chrom"}, \code{"exon_strand"},
    \code{"cds_id"}, \code{"cds_name"},
    \code{"cds_chrom"}, \code{"cds_strand"} and \code{"exon_rank"}.
  }
  \item{columns}{
    Columns to include in the output.
    Must be \code{NULL} or a character vector with values in the
    above list of valid names. With the following restrictions:
    \itemize{
    \item \code{"tx_chrom"} and \code{"tx_strand"} are not allowed
          for \code{transcripts}.
    \item \code{"exon_chrom"} and \code{"exon_strand"} are not allowed
          for \code{exons}.
    \item \code{"cds_chrom"} and \code{"cds_strand"} are not allowed
          for \code{cds}.
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
  exons(txdb, vals=list(exon_id=1), columns=c("exon_id", "tx_name"))
  exons(txdb, vals=list(tx_name="uc009vip.1"), columns=c("exon_id", "tx_name"))
}
