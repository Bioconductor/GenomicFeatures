\name{makeTranscriptDbFromUCSC}

\alias{supportedUCSCtables}
\alias{getChromInfoFromUCSC}
\alias{makeTranscriptDbFromUCSC}

\title{
  Make a TxDb object from annotations available at the
  UCSC Genome Browser
}
\description{
  The \code{makeTranscriptDbFromUCSC} function allows the user
  to make a \link{TxDb} object from transcript annotations
  available at the UCSC Genome Browser.
}
\usage{
supportedUCSCtables()

getChromInfoFromUCSC(
        genome,
        goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")

makeTranscriptDbFromUCSC(
        genome="hg19",
        tablename="knownGene",
        transcript_ids=NULL,
        circ_seqs=DEFAULT_CIRC_SEQS,
        url="http://genome.ucsc.edu/cgi-bin/",
        goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
        miRBaseBuild=NA)
}
\arguments{
  \item{genome}{genome abbreviation used by UCSC and obtained by
    \code{\link[rtracklayer]{ucscGenomes}()[ , "db"]}.
    For example: \code{"hg19"}.}
  \item{tablename}{name of the UCSC table containing the transcript
    annotations to retrieve. Use the \code{supportedUCSCtables} utility
    function to get the list of supported tables. Note that not all
    tables are available for all genomes.}
  \item{transcript_ids}{optionally, only retrieve transcript
    annotation data for the specified set of transcript ids.
    If this is used, then the meta information displayed for the
    resulting \link{TxDb} object will say 'Full dataset: no'.
    Otherwise it will say 'Full dataset: yes'.}
  \item{circ_seqs}{a character vector to list out which chromosomes
    should be marked as circular.}
  \item{url,goldenPath_url}{use to specify the location of an
    alternate UCSC Genome Browser.}
  \item{miRBaseBuild}{specify the string for the appropriate build
    Information from mirbase.db to use for microRNAs.  This can be
    learned by calling \code{supportedMiRBaseBuildValues}.  By default,
    this value will be set to \code{NA}, which will inactivate the
    \code{microRNAs} accessor.}
}
\details{
  \code{makeTranscriptDbFromUCSC} is a convenience function that feeds
  data from the UCSC source to the lower level \code{\link{makeTranscriptDb}}
  function.
  See \code{?\link{makeTranscriptDbFromBiomart}} for a similar function
  that feeds data from a BioMart database.
}

\value{A \link{TxDb} object.}

\author{
  M. Carlson and H. Pages
}

\seealso{
  \code{\link[rtracklayer]{ucscGenomes}},
  \code{\link{DEFAULT_CIRC_SEQS}},
  \code{\link{makeTranscriptDbFromBiomart}},
  \code{\link{makeTranscriptDbFromGFF}},
  \code{\link{makeTranscriptDb}},
  \code{\link{supportedMiRBaseBuildValues}}
}

\examples{
## Display the list of genomes available at UCSC:
library(rtracklayer)
ucscGenomes()[ , "db"]

## Display the list of tables supported by makeTranscriptDbFromUCSC():
supportedUCSCtables()

\dontrun{
## Retrieving a full transcript dataset for Yeast from UCSC:
txdb1 <- makeTranscriptDbFromUCSC(genome="sacCer2", tablename="ensGene")
}

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
