\name{nearest-methods}

\alias{nearest-methods}

\alias{distance,GenomicRanges,TxDb-method}

\title{Finding the nearest genomic range neighbor in a TxDb}

\description{
  The \code{distance} methods for TxDb objects and subclasses.
}

\usage{
\S4method{distance}{GenomicRanges,TxDb}(x, y, ignore.strand=FALSE,
    ..., id, type=c("gene", "tx", "exon", "cds"))
}

\arguments{
  \item{x}{The query \link{GenomicRanges} instance.
  }
  \item{y}{For \code{distance}, a \link{TxDb} instance. The
   \code{id} is used to extract ranges from the \link{TxDb}
   which are then used to compute the distance from \code{x}.
  }
  \item{id}{A \code{character} vector the same length as \code{x}.
    The \code{id} must be identifiers in the \link{TxDb} object.
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
       extracted from the \link{TxDb} object \code{y}. Values in 
       \code{id} are matched to one of \sQuote{gene_id}, \sQuote{tx_id}, 
       \sQuote{exon_id} or \sQuote{cds_id} identifiers in the \link{TxDb}
       and the corresponding ranges are extracted. The \code{type} argument 
       specifies which identifier is represented in \code{id}. The extracted 
       ranges are used in the distance calculation with the ranges in 
       \code{x}.

       The method returns \code{NA} values when the genomic region defined
       by \code{id} cannot be collapsed into a single range (e.g.,
       when a gene spans multiple chromosomes) or if the \code{id}
       is not found in \code{y}.

       The behavior of \code{distance} with respect to zero-width ranges
       has changed in Bioconductor 2.12. See the man page \code{?distance} 
       in IRanges for details.
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
  distance(gr, txdb, id=c("10000", "23000"), type="cds")

  ## The id's must be in the appropriate order with respect to 'x'.
  distance(gr, txdb, id=c("4", "4097"), type="tx")

  ## 'id' "4" is on chr2L and "4097" is on chr2R. 
  transcripts(txdb, filter=list(tx_id=c("4", "4097")))

  ## If we reverse the 'id' the chromosomes are incompatable with gr.
  distance(gr, txdb, id=c("4097", "4"), type="tx")

  ## distance() compares each 'x' to the corresponding 'y'.
  ## If an 'id' is not found in the TxDb 'y' will not 
  ## be the same lenth as 'x' and an error is thrown.
  \dontrun{ 
  distance(gr, txdb, id=c("FBgn0000008", "INVALID"), type="gene") ## will fail
  } 
}

\keyword{utilities}
