\name{microRNAs}

\alias{microRNAs}
\alias{microRNAs,TxDb-method}
\alias{tRNAs}
\alias{tRNAs,TxDb-method}

\title{
  Extract microRNA or tRNA genomic ranges from an object
}

\description{
  Generic functions to extract microRNA or tRNA genomic ranges from an object.
  This page documents the methods for \link{TxDb} objects only.
}

\usage{
microRNAs(x)
\S4method{microRNAs}{TxDb}(x)

tRNAs(x)
\S4method{tRNAs}{TxDb}(x)
}

\arguments{ 
  \item{x}{
    A \link{TxDb} object.
  }
}

\value{
  A \link[GenomicRanges]{GRanges} object.
}

\author{
  M. Carlson
}

\seealso{
  \itemize{
    \item \code{\link{transcripts}}, \code{\link{transcriptsBy}}, and
          \code{\link{transcriptsByOverlaps}} for the core genomic features
          extractors.
    \item The \link{TxDb} class.
  }
}

\examples{
\dontrun{library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(mirbase.db)
microRNAs(TxDb.Hsapiens.UCSC.hg19.knownGene)
}
}

\keyword{methods}
