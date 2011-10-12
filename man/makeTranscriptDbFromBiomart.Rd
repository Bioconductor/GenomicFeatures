\name{makeTranscriptDbFromBiomart}

\alias{makeTranscriptDbFromBiomart}
\alias{getChromInfoFromBiomart}

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
getChromInfoFromBiomart(biomart="ensembl",
                        dataset="hsapiens_gene_ensembl")

makeTranscriptDbFromBiomart(biomart="ensembl",
                            dataset="hsapiens_gene_ensembl",
                            transcript_ids=NULL,
                            circ_seqs=DEFAULT_CIRC_SEQS)
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
  \item{circ_seqs}{a character vector to list out which chromosomes
    should be marked as circular.}
}

\details{
  \code{makeTranscriptDbFromBiomart} is a convenience function that feeds
  data from a BioMart database to the lower level \code{\link{makeTranscriptDb}}
  function.
  See \code{?\link{makeTranscriptDbFromUCSC}} for a similar function
  that feeds data from the UCSC source.

  BioMart databases that are known to have compatible transcript annotations are:
  \itemize{
    \item the most recent ensembl: ENSEMBL GENES (SANGER UK)
    \item the most recent bacterial_mart: ENSEMBL BACTERIA (EBI UK)
    \item the most recent fungal_mart: ENSEMBL FUNGAL (EBI UK)
    \item the most recent metazoa_mart: ENSEMBL METAZOA (EBI UK)
    \item the most recent plant_mart: ENSEMBL PLANT (EBI UK)
    \item the most recent protist_mart: ENSEMBL PROTISTS (EBI UK)
    \item the most recent ensembl_expressionmart: EURATMART (EBI UK)
  }
  Not all annotations will have CDS information.
}

\value{A \link{TranscriptDb} object.}

\author{
  M. Carlson and H. Pages
}

\seealso{
  \code{\link[biomaRt]{listMarts}},
  \code{\link[biomaRt]{useMart}},
  \code{\link[biomaRt]{listDatasets}},
  \code{\link{DEFAULT_CIRC_SEQS}},
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
    "ENST00000268655",
    "ENST00000313243",
    "ENST00000341724",
    "ENST00000400839",
    "ENST00000400840",
    "ENST00000435657",
    "ENST00000478783"
)
txdb <- makeTranscriptDbFromBiomart(transcript_ids=transcript_ids)
txdb  # note that these annotations match the GRCh37 genome assembly
}
