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
\alias{seqinfo,TranscriptDb-method}
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
      \code{seqinfo(x)}:
      Gets the information about the underlying sequences as a
      \link[GenomicRanges]{Seqinfo} object.
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
  \link[GenomicRanges]{Seqinfo-class},
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

seqinfo(txdb)
seqlevels(txdb)  # shortcut for 'seqlevels(seqinfo(txdb))'
seqlengths(txdb)  # shortcut for 'seqlengths(seqinfo(txdb))'
isCircular(txdb)  # shortcut for 'isCircular(seqinfo(txdb))'
names(which(isCircular(txdb)))

txdump <- as.list(txdb)
txdump
txdb1 <- do.call(makeTranscriptDb, txdump)
stopifnot(identical(as.list(txdb1), txdump))
}
