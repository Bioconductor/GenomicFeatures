\name{features}

\alias{features}
\alias{features,FeatureDb-method}

\title{
  Extract simple features from a FeatureDb object
}

\description{
  WARNING: The FeatureDb/makeFeatureDbFromUCSC/features code base is
  no longer actively maintained and FeatureDb-related functionalities
  might get deprecated in the near future. Please use
  \code{\link{makeFeatureDbFromUCSC}} for a convenient way to
  import transcript annotations from UCSC online resources into
  Bioconductor.

  Generic function to extract genomic features from a FeatureDb object.
}

\usage{
features(x)
\S4method{features}{FeatureDb}(x)
}

\arguments{  
  \item{x}{
    A \link{FeatureDb} object.
  }
}


\value{ a GRanges object }

\author{
  M. Carlson
}

\seealso{
  \link{FeatureDb}
}

\examples{
  fdb <- loadDb(system.file("extdata", "FeatureDb.sqlite",
                                   package="GenomicFeatures"))
  features(fdb)
}
