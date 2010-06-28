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
  \item{biomart}{which BioMart database to use.
    Get the list of all available BioMart databases with the
    \code{\link[biomaRt]{listMarts}} function from the biomaRt
    package. See the details section below for a list of BioMart
    databases with compatible transcript annotations.}
  \item{dataset}{which dataset from BioMart. For example:
    \code{"hsapiens_gene_ensembl"}, \code{"mmusculus_gene_ensembl"},
    \code{"dmelanogaster_gene_ensembl"}, \code{"celegans_gene_ensembl"},
    \code{"scerevisiae_gene_ensembl"}, etc in the ensembl database.
    See the examples section below for how to discover which datasets
    are available in a given BioMart database.}
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

  As of June 28, 2010, the BioMart databases with compatible transcript
  annotations are:
  \itemize{
    \item ensembl: ENSEMBL GENES 58 (SANGER UK)
    \item bacterial_mart_5: ENSEMBL BACTERIA 5 (EBI UK)
    \item fungal_mart_5: ENSEMBL FUNGAL 5 (EBI UK)
    \item metazoa_mart_5: ENSEMBL METAZOA 5 (EBI UK)
    \item plant_mart_5: ENSEMBL PLANT 5 (EBI UK)
    \item protist_mart_5: ENSEMBL PROTISTS 5 (EBI UK)
    \item ensembl_expressionmart_48: EURATMART (EBI UK)
    \item Ensembl56: PANCREATIC EXPRESSION DATABASE (INSTITUTE OF CANCER UK)
    \item ENSEMBL_MART_ENSEMBL: GRAMENE 30 ENSEMBL GENES (CSHL/CORNELL US)
  }
  Only ensembl and Ensembl56 have CDS information.
}

\value{A \link{TranscriptDb} object.}

\author{
  M. Carlson and H. Pages
}

\seealso{
  \code{\link[biomaRt]{listMarts}},
  \code{\link[biomaRt]{useMart}},
  \code{\link[biomaRt]{listDatasets}},
  \code{\link{makeTranscriptDbFromUCSC}},
  \code{\link{makeTranscriptDb}}
}

\examples{
## Discover which datasets are available in the "ensembl" and
## "plant_mart_5" BioMart databases:
library(biomaRt)
listDatasets(useMart("ensembl"))
listDatasets(useMart("plant_mart_5"))

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
