\name{FeatureDb-class}

\alias{FeatureDb-class}
\alias{class:FeatureDb}
\alias{FeatureDb}

\title{FeatureDb objects}

\description{
  The FeatureDb class is a generic container for storing
  genomic locations of an arbitrary type of genomic features.
  
  See \code{?\link{TxDb}} for a container for storing transcript
  annotations.

  See \code{?\link{makeFeatureDbFromUCSC}} for a convenient way to
  make FeatureDb objects from BioMart online resources.
}

\section{Methods}{
  In the code snippets below, \code{x} is a FeatureDb object.

  \describe{
    \item{}{
      \code{metadata(x)}:
      Return \code{x}'s metadata in a data frame.
    }
  }
}

\author{Marc Carlson}

\seealso{
  \itemize{
    \item The \link{TxDb} class for storing transcript annotations.
    \item \code{\link{makeFeatureDbFromUCSC}} for a convenient way to
          make a FeatureDb object from UCSC online resources.
    \item \code{\link{saveDb}} and \code{\link{loadDb}} for
          saving and loading the database content of a FeatureDb object.
    \item \code{\link{features}} for how to extract genomic features
          from a FeatureDb object.
  }
}

\examples{
fdb_file <- system.file("extdata", "FeatureDb.sqlite",
                        package="GenomicFeatures")
fdb <- loadDb(fdb_file)
fdb
}

\keyword{methods}
\keyword{classes}
