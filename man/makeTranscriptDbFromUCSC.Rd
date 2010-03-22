\name{makeTranscriptDbFromUCSC}

\alias{supportedUCSCtables}
\alias{makeTranscriptDbFromUCSC}

\title{
  Making a TranscriptDb object from annotations available at the
  UCSC genome browser
}
\description{
  The \code{makeTranscriptDbFromUCSC} function allows the user
  to make a \link{TranscriptDb} object from transcript annotations
  available at the UCSC genome browser.
}
\usage{
  supportedUCSCtables()
  makeTranscriptDbFromUCSC(genome="hg18",
                           tablename="knownGene",
                           transcript_ids=NULL)
}
\arguments{
  \item{genome}{genome abbreviation used by UCSC and obtained by
    \code{\link[rtracklayer]{ucscGenomes}()[ , "db"]}.
    For example: \code{"hg18"}.}
  \item{tablename}{name of the UCSC table containing the transcript
    annotations to retrieve. Use the \code{supportedUCSCtables} utility
    function to get the list of supported tables. Note that not all
    tables are available for all genomes.}
  \item{transcript_ids}{optionally, only retrieve transcript
    annotation data for the specified set of transcript ids.
    If this is used, then the meta information displayed for the
    resulting \link{TranscriptDb} object will say 'Full dataset: no'.
    Otherwise it will say 'Full dataset: yes'.}
}
\details{
  \code{makeTranscriptDbFromUCSC} is a convenience function that feeds
  data from the UCSC source to the lower level \code{\link{makeTranscriptDb}}
  function.
  See \code{?\link{makeTranscriptDbFromBiomart}} for a similar function
  that feeds data from a BioMart database.
}

\value{A \link{TranscriptDb} object.}

\author{
  M. Carlson and H. Pages
}

\seealso{
  \code{\link[rtracklayer]{ucscGenomes}},
  \code{\link{makeTranscriptDbFromBiomart}},
  \code{\link{makeTranscriptDb}}
}

\examples{
## Display the list of genomes available at UCSC:
library(rtracklayer)
ucscGenomes()[ , "db"]

## Display the list of tables supported by makeTranscriptDbFromUCSC():
supportedUCSCtables()

## Retrieving a full transcript dataset for Yeast from UCSC:
txdb1 <- makeTranscriptDbFromUCSC(genome="sacCer2", tablename="ensGene")
txdb1

## Retrieving an incomplete transcript dataset for Mouse from UCSC
## (only transcripts linked to Entrez Gene ID 22290):
transcript_ids <- c(
    "uc009uzf.1",
    "uc009uzg.1",
    "uc009uzh.1",
    "uc009uzi.1",
    "uc009uzj.1"
)

txdb2 <- makeTranscriptDbFromUCSC(genome="mm9", tablename="knownGene",
                                  transcript_ids=transcript_ids)
txdb2
}
