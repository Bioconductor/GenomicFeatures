\name{makeFeatureDbFromUCSC}

\alias{supportedUCSCFeatureDbTracks}
\alias{makeFeatureDbFromUCSC}

\title{
  Making a FeatureDb object from annotations available at the
  UCSC Genome Browser
}
\description{
  The \code{makeFeatureDbFromUCSC} function allows the user
  to make a \link{FeatureDb} object from simple annotation tracks at
  UCSC.  The tracks in question must (at a minimum) have a start, end and
  a chromosome affiliation in order to be made into a \link{FeatureDb}.
}
\usage{
supportedUCSCFeatureDbTracks()

makeFeatureDbFromUCSC(
        genome="hg18",
        tablename="oreganno",
        columns = c(id="character",name="character"),
        url="http://genome.ucsc.edu/cgi-bin/",
        goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
}
\arguments{
  \item{genome}{genome abbreviation used by UCSC and obtained by
    \code{\link[rtracklayer]{ucscGenomes}()[ , "db"]}.
    For example: \code{"hg18"}.}
  \item{tablename}{name of the UCSC table containing the transcript
    annotations to retrieve. Use the \code{supportedUCSCtables} utility
    function to get the list of supported tables. Note that not all
    tables are available for all genomes.}
  \item{columns}{a named character vector to list out the names and
        types of the other columns that the downloaded track should have.}
  \item{url,goldenPath_url}{use to specify the location of an
    alternate UCSC Genome Browser.}
}
\details{
  \code{makeFeatureDbFromUCSC} is a convenience function that feeds
  data from the UCSC source to a lower level function.
  % \code{\link{makeFeatureDb}}
  % function.
  % See \code{?\link{makeFeatureDbFromBiomart}} for a similar function
  % that feeds data from a BioMart database.
}

\value{A \link{FeatureDb} object.}

\author{
  M. Carlson and H. Pages
}

\seealso{
  \code{\link[rtracklayer]{ucscGenomes}},
  % \code{\link{makeFeatureDbFromBiomart}},
  % \code{\link{makeFeatureDb}}
}

\examples{
## Display the list of genomes available at UCSC:
library(rtracklayer)
ucscGenomes()[ , "db"]

## Display the list of tables supported by makeFeatureDbFromUCSC():
supportedUCSCFeatureDbTracks()

## Retrieving a full transcript dataset for Yeast from UCSC:
txdb <- makeFeatureDbFromUCSC(genome="mm9", tablename="oreganno")
txdb
}
