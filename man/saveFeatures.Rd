\name{saveFeatures}
\alias{saveFeatures}
\alias{loadFeatures}
\title{
  Methods to save and load the database contents for a Transcript Object.
}
\description{
  These methods provide a way to dump a TranscriptDb object to
  an SQLite file, and to recreate that object the saved file.
}
\usage{
  saveFeatures(x, file)
  loadFeatures(file)
}
\arguments{
  \item{x}{a transcripts object, which contains a connection to a DB.}
  \item{file}{A SQLite Database filename.}
}
\value{
  For \code{loadFeatures} only, a \link{TranscriptDb}
  object is returned.
}
\author{
  M. Carlson
}
\seealso{
  \link{TranscriptDb}
}
\examples{
  txdb <- 
    loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                             package = "GenomicFeatures"))
  txdb
}
