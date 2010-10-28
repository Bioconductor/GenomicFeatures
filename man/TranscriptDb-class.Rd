\name{TranscriptDb-class}

\alias{GenomicFeatures-class}
\alias{class:GenomicFeatures}
\alias{GenomicFeatures}

\alias{TranscriptDb-class}
\alias{class:TranscriptDb}
\alias{TranscriptDb}
\alias{FeatureDb-class}
\alias{class:FeatureDb}
\alias{FeatureDb}


\alias{metadata,TranscriptDb-method}
\alias{metadata,FeatureDb-method}
\alias{seqnames,TranscriptDb-method}
\alias{seqlengths,TranscriptDb-method}
\alias{isCircular,TranscriptDb-method}
\alias{show,TranscriptDb-method}
\alias{show,FeatureDb-method}
\alias{as.list,TranscriptDb-method}

\title{TranscriptDb objects}

\description{
  The  TranscriptDb class is a container for storing transcript
  annotations.
  The FeatureDb class is a container for storing more generic
  GenomicFeature annotations.

  See \code{?\link{makeTranscriptDbFromUCSC}} and
  \code{?\link{makeTranscriptDbFromBiomart}} for making a
  TranscriptDb object from the UCSC or BioMart sources.

  See \code{?\link{makeFeatureDbFromUCSC}} for making a
  FeatureDb object from the UCSC or BioMart sources.
  
  See \code{?\link{saveFeatures}} and \code{?\link{loadFeatures}} for
  saving and loading the database contents of a TranscriptDb or FeatureDb
  object.
}

\section{Methods}{
  In the code snippets below, \code{x} is a TranscriptDb object. For the
  metadata and show methods, there is also support for FeatureDb objects.

  \describe{
    \item{}{
      \code{metadata(x)}:
      Returns \code{x}'s metadata in a data frame.
    }
    \item{}{
      \code{seqnames(x)}:
      Returns the names of all chromosomes in a character vector.
      Note that "all" here means at least the chromosomes that have
      features. But some TranscriptDb objects (in particular those
      created with \code{\link{makeTranscriptDbFromUCSC}})
      are storing the names of all the chromosomes forming the genome.
    }
    \item{}{
      \code{seqlengths(x)}:
      Returns the lengths of the chromosomes \emph{or NAs}
      in an integer vector named with \code{seqnames(x)}.
      Either \emph{all} elements are NAs or \emph{none} is,
      depending on the avaibality of this information at creation time.
    }
    \item{}{
      \code{isCircular(x)}:
      Returns the circularity flag of the chromosomes \emph{or NAs}
      in a logical vector named with \code{seqnames(x)}.
      Either \emph{all} elements are NAs or \emph{none} is,
      depending on the avaibality of this information at creation time.
    }
    \item{}{
      \code{as.list(x)}:
      Dumps the entire db into a list of data frames \code{txdump} that can
      be used in \code{do.call(makeTranscriptDb, txdump)} to make the db
      again with no loss of information.
      Note that the transcripts are dumped in the same order in all the
      data frames.
    }
  }

  See \code{?\link{transcripts}}, \code{?\link{transcriptsByOverlaps}},
  \code{?\link{id2name}} and \code{?\link{transcriptsBy}} for other useful
  operations on TranscriptDb objects.
}

\author{H. Pages}

\seealso{
  \code{\link{makeTranscriptDbFromUCSC}},
  \code{\link{makeTranscriptDbFromBiomart}},
  \code{\link{loadFeatures}},
  \code{\link{transcripts}},
  \code{\link{transcriptsByOverlaps}},
  \code{\link{id2name}}, \code{\link{transcriptsBy}}
}

\examples{
txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                         package="GenomicFeatures")
txdb <- loadFeatures(txdb_file)
txdb

seqnames(txdb)
seqlengths(txdb)
names(which(isCircular(txdb)))

txdump <- as.list(txdb)
txdump
txdb1 <- do.call(makeTranscriptDb, txdump)
stopifnot(identical(as.list(txdb1), txdump))
}
