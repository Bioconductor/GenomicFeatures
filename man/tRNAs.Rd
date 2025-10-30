\name{tRNAs}

\alias{tRNAs}
\alias{tRNAs,TxDb-method}

\title{
  Extract tRNA genomic ranges from an object
}

\description{
  WARNING: The code base for \code{tRNAs()} is no longer actively
  maintained and the function might get deprecated in the near future.

  The \code{tRNAs()} function extracts tRNA genomic ranges from a
  \link{TxDb} object.
}

\usage{
tRNAs(x)
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

\keyword{methods}
