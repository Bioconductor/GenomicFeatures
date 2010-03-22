\name{makeTranscriptDbFromBiomart}

\alias{makeTranscriptDbFromBiomart}

\title{
  Making a TranscriptDb object from annotations available on a
  BioMart database
}
\description{
  The \code{makeTranscriptDbFromBiomart} function allows the user
  to make a \link{TranscriptDb} object from transcript annotations
  available on a BioMart database.
}
\usage{
  makeTranscriptDbFromBiomart(biomart="ensembl",
                              dataset="hsapiens_gene_ensembl",
                              transcript_ids=NULL)
}
\arguments{
  \item{biomart}{which BioMart database to use. Only \code{"ensembl"}
    (the default) is supported for now.}
  \item{dataset}{which dataset from BioMart. For example:
    \code{"hsapiens_gene_ensembl"}, \code{"mmusculus_gene_ensembl"},
    \code{"dmelanogaster_gene_ensembl"}, \code{"celegans_gene_ensembl"},
    \code{"scerevisiae_gene_ensembl"}, etc. See the examples section
    below for how to discover which datasets are available in a given
    BioMart database.}
  \item{transcript_ids}{optionally, only retrieve transcript
    annotation data for the specified set of transcript ids.
    If this is used, then the meta information displayed for the
    resulting \link{TranscriptDb} object will say 'Full dataset: no'.
    Otherwise it will say 'Full dataset: yes'.}
}
\details{
  \code{makeTranscriptDbFromBiomart} is a convenience function that feeds
  data from a BioMart database to the lower level \code{\link{makeTranscriptDb}}
  function.
  See \code{?\link{makeTranscriptDbFromUCSC}} for a similar function
  that feeds data from the UCSC source.
}

\value{A \link{TranscriptDb} object.}

\author{
  M. Carlson and H. Pages
}

\seealso{
  \code{\link[biomaRt]{useMart}},
  \code{\link[biomaRt]{listDatasets}},
  \code{\link{makeTranscriptDbFromUCSC}},
  \code{\link{makeTranscriptDb}}
}

\examples{
## Discover which datasets are available in the "ensembl" BioMart
## database:
library(biomaRt)
listDatasets(useMart("ensembl"))

## Retrieving an incomplete transcript dataset for Human from the
## "ensembl" BioMart database:
transcript_ids <- c(
    "ENST00000400839",
    "ENST00000400840",
    "ENST00000478783",
    "ENST00000435657",
    "ENST00000268655",
    "ENST00000313243",
    "ENST00000341724"
)
txdb <- makeTranscriptDbFromBiomart(transcript_ids=transcript_ids)
txdb  # note that these annotations match the GRCh37 genome assembly
}
