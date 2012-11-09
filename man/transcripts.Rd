\name{transcripts}

\alias{transcripts}
\alias{transcripts,data.frame-method}
\alias{transcripts,TranscriptDb-method}
\alias{exons}
\alias{exons,data.frame-method}
\alias{exons,TranscriptDb-method}
\alias{cds}
\alias{cds,TranscriptDb-method}
\alias{promoters,TranscriptDb-method}
\alias{microRNAs}
\alias{microRNAs,TranscriptDb-method}
\alias{tRNAs}
\alias{tRNAs,TranscriptDb-method}

\title{
  Extract genomic features from an object
}

\description{
  Generic functions to extract genomic features from an object.
  This page documents the methods for \link{TranscriptDb} objects only.
}

\usage{
transcripts(x, ...)
\S4method{transcripts}{TranscriptDb}(x, vals=NULL, columns=c("tx_id", "tx_name"))

exons(x, ...)
\S4method{exons}{TranscriptDb}(x, vals=NULL, columns="exon_id")

cds(x, ...)
\S4method{cds}{TranscriptDb}(x, vals=NULL, columns="cds_id")

#promoters(x, upstream=2000, downstream=200, ...)
\S4method{promoters}{TranscriptDb}(x, upstream=2000, downstream=200, ...)

microRNAs(x)
\S4method{microRNAs}{TranscriptDb}(x)

tRNAs(x)
\S4method{tRNAs}{TranscriptDb}(x)
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
    If the vector is named, those names are used for the corresponding
    column in the element metadata of the returned object. 
  }
  \item{upstream}{
    For \code{promoters()} : An \code{integer(1)} value indicating the
    number of bases upstream from the transcription start site. For
    additional details see ?\code{'promoters,GRanges-method'}.
  }
  \item{downstream}{
    For \code{promoters()} : An \code{integer(1)} value indicating the 
    number of bases downstream from the transcription start site. For 
    additional details see ?\code{'promoters,GRanges-method'}.
  }
}

\details{
  These are the main functions for extracting transcript information
  from a \link{TranscriptDb} object. With the exception of
  \code{microRNAs}, these methods can restrict the output based on
  categorical information. To restrict the output based on interval
  information, use the \code{\link{transcriptsByOverlaps}},
  \code{\link{exonsByOverlaps}}, and \code{\link{cdsByOverlaps}}
  functions.

  The \code{promoters()} function computes user-defined promoter regions
  for the transcripts in a \code{TranscriptDb} object. The return object is a 
  \code{GRanges} of promoter regions around the transcription start 
  site the span of which is defined by \code{upstream} and \code{downstream}.
  For additional details on how the promoter range is computed and the
  handling of \code{+} and \code{-} strands see 
  ?\code{'promoters,GRanges-method'}.
}

\value{ a GRanges object }

\author{
  M. Carlson, P. Aboyoun and H. Pages
}

\seealso{
  \itemize{
    \item \code{\link{transcriptsBy}} and \code{\link{transcriptsByOverlaps}}
          for more ways to extract genomic features
          from a \link{TranscriptDb} object.
    \item \link[GenomicFeatures]{select-methods} for how to use the
          simple "select" interface to extract information from a
          \link{TranscriptDb} object.
    \item \code{\link{id2name}} for mapping \link{TranscriptDb} internal ids
          to external names for a given feature type.
    \item The \link{TranscriptDb} class.
  }
}

\examples{
## transcripts() and exons() :
  txdb <- loadDb(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                   package="GenomicFeatures"))
  vals <- list(tx_chrom = c("chr3", "chr5"), tx_strand = "+")
  transcripts(txdb, vals)
  exons(txdb, vals=list(exon_id=1), columns=c("exon_id", "tx_name"))
  exons(txdb, vals=list(tx_name="uc009vip.1"), columns=c("exon_id",
      "tx_name"))

## microRNAs() :
\dontrun{library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(mirbase.db)
microRNAs(TxDb.Hsapiens.UCSC.hg19.knownGene)
}

## promoters() :
head(promoters(txdb, 100, 50))
}

\keyword{methods}
