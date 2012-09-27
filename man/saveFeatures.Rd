\name{saveFeatures}

\alias{saveFeatures}
\alias{saveFeatures,TranscriptDb-method}
\alias{saveFeatures,FeatureDb-method}
\alias{loadFeatures}
\title{
  Methods to save and load the database contents for a TranscriptDb
  or FeatureDb object.
}
\description{
  These methods provide a way to dump a \link{TranscriptDb} or
  \link{FeatureDb} object to an SQLite file, and to recreate that object
  from the saved file.

  However, these methods are now deprecated and have been replaced by
  \code{\link[AnnotationDbi]{saveDb}} and \code{\link[AnnotationDbi]{loadDb}}.

  Users are encouraged to switch to those other methods as the methods
  documented here will soon be defunct.
}
\usage{
  saveFeatures(x, file)
  loadFeatures(file)
}
\arguments{
  \item{x}{A \link{TranscriptDb} or \link{FeatureDb} object.}
  \item{file}{An SQLite Database filename.}
}
\value{
  For \code{loadFeatures} only, a \link{TranscriptDb} or \link{FeatureDb}
  object is returned.
}
\author{
  M. Carlson
}
\seealso{
  \code{\link[AnnotationDbi]{saveDb}},
  \link{TranscriptDb},
  \link{FeatureDb}
}
\examples{
\dontrun{
  txdb <- 
    loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                             package = "GenomicFeatures"))
  txdb
}
}
