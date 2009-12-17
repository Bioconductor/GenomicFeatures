\name{saveFeatures}
\alias{saveFeatures}
\alias{loadFeatures}
\title{
  Methods to save and load the database contents for a Transcript Object.
}
\description{
  These methods provide a way to save the contents of a generated
  Transcripts Object Database for for future use (as a SQLIte database),
  and to recreate that object from a source Database file.
}
\usage{
  saveFeatures(x, file)
  loadFeatures(file)
}
\arguments{
  \item{x}{a transcripts object, which contains a connection to a DB.}
  \item{file}{A SQLite Database filename.}
}
\value{For loadFeatures only, a Transcripts object is returned.}
\author{
  M. Carlson
}

\examples{
\dontrun{
  ##Have a Transcript Object called tx?
  ##You can save it:
  saveFeatures(tx, file="HG18.sqlite")
  ##Then later you can quickly  re-constitute it like this:
  tx <- loadFeatures("HG18.sqlite")
}
}