\name{TranscriptDb-class}

\alias{GenomicFeatures-class}
\alias{class:GenomicFeatures}
\alias{GenomicFeatures}

\alias{TranscriptDb-class}
\alias{class:TranscriptDb}
\alias{TranscriptDb}

\alias{show,TranscriptDb-method}
\alias{as.list,TranscriptDb-method}

\title{TranscriptDb objects}

\description{
  The TranscriptDb class is a container for storing transcript annotations.

  See \code{?\link{makeTranscriptDbFromUCSC}} and
  \code{?\link{makeTranscriptDbFromBiomart}} for making a
  TranscriptDb object from the UCSC or BioMart sources.

  See \code{?\link{saveFeatures}} and \code{?\link{loadFeatures}} for
  saving and loading the database contents of a TranscriptDb object.
}

\section{Methods}{
  In the code snippets below, \code{x} is a TranscriptDb object.

  \describe{
    \item{}{
      \code{as.list(x)}:
      Dumps the entire db into a list of data frames \code{txdump} that can
      be used in \code{do.call(makeTranscriptDb, txdump)} to make the db
      again with no loss of information.
      Note that the transcripts are dumped in the same order in all the
      data frames.
    }
  }

  See \code{?\link{transcripts}}, \code{?\link{transcriptsByRanges}}
  and \code{?\link{bindTranscripts}} for other useful operations on
  TranscriptDb objects.
}

\author{H. Pages}

\seealso{
  \code{\link{makeTranscriptDbFromUCSC}},
  \code{\link{makeTranscriptDbFromBiomart}},
  \code{\link{loadFeatures}},
  \code{\link{transcripts}},
  \code{\link{transcriptsByRanges}},
  \code{\link{bindTranscripts}}
}

\examples{
txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                         package="GenomicFeatures")
txdb <- loadFeatures(txdb_file)
txdb

txdump <- as.list(txdb)
txdump
txdb1 <- do.call(makeTranscriptDb, txdump)
stopifnot(identical(as.list(txdb1), txdump))
}
