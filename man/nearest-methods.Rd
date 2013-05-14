\name{nearest-methods}

\alias{nearest-methods}

\alias{distance,GenomicRanges,TranscriptDb-method}

\title{Finding the nearest genomic range neighbor in a TranscriptDb}

\description{
  The \code{distance} methods for TranscriptDb objects and subclasses.
}

\usage{
\S4method{distance}{GenomicRanges,TranscriptDb}(x, y, ignore.strand=FALSE,
    ..., id, type=c("gene", "tx", "exon", "cds"))
}

\arguments{
  \item{x}{The query \link{GenomicRanges} instance.
  }
  \item{y}{For \code{distance}, a \link{TranscriptDb} instance. The
   \code{id} is used to extract ranges from the \link{TranscriptDb}
   which are then used to compute the distance from \code{x}.
  }
  \item{id}{A \code{character} vector the same length as \code{x}.
    The \code{id} must be identifiers in the TranscriptDb object.
    \code{type} indicates what type of identifier \code{id} is. 
  }
  \item{type}{A \code{character(1)} describing the \code{id}. 
    Must be one of \sQuote{gene}, \sQuote{tx}, \sQuote{exon} or
    \sQuote{cds}.
  }
  \item{ignore.strand}{
    A \code{logical} indicating if the strand of the ranges
    should be ignored. When \code{TRUE}, strand is set to \code{'+'}.
  }
  \item{...}{Additional arguments for methods.
  }
}

\details{
  \itemize{
    \item{distance: }{
       Returns the distance for each range in \code{x} to the range 
       extracted from the TranscriptDb object \code{y}. Values in 
       \code{id} are matched to one of \sQuote{gene_id}, \sQuote{tx_id}, 
       \sQuote{exon_id} or \sQuote{cds_id} identifiers in the TranscriptDb 
       and the corresponding ranges are extracted. The \code{type} argument 
       specifies which identifier is represented in \code{id}. The extracted 
       ranges are used in the distance calculation with the ranges in 
       \code{x}.

       The behavior of \code{distance} has changed in Bioconductor 2.12.
       See the man page \code{?distance} in IRanges for details.
    }
  }
}

\value{
  For \code{distance}, an integer vector of distances between the ranges
  in \code{x} and \code{y}.
}

\author{Valerie Obenchain <vobencha@fhcrc.org>}


\seealso{
  \itemize{
    \item \link[IRanges]{nearest-methods} man page in IRanges.
    \item \link[GenomicRanges]{nearest-methods} man page in GenomicRanges.
  }
}

\examples{
  ## -----------------------------------------------------------
  ## distance()
  ## -----------------------------------------------------------

  library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
  txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
  gr <- GRanges(c("chr2L", "chr2R"),
                IRanges(c(100000, 200000),  width=100))
  distance(gr, txdb, id=c("FBgn0259717", "FBgn0261501"), type="gene")
  distance(gr, txdb, id=c("4", "4097"), type="tx")
  distance(gr, txdb, id=c("10000", "23000"), type="cds")

  ## distance() compares each 'x' to the corresponding 'y'.
  ## If an 'id' is not found in the TranscriptDb 'y' will not 
  ## be the same lenth as 'x' and an error is thrown.
  \dontrun{ 
  distance(gr, txdb, id=c("FBgn0000008", "INVALID"), type="gene") ## will fail
  } 
}

\keyword{utilities}
