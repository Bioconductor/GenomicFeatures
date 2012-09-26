\name{TranscriptDb-class}

\alias{TranscriptDb-class}
\alias{class:TranscriptDb}
\alias{TranscriptDb}
\alias{FeatureDb-class}
\alias{class:FeatureDb}
\alias{FeatureDb}

\alias{metadata,TranscriptDb-method}
\alias{metadata,FeatureDb-method}
\alias{seqinfo,TranscriptDb-method}
\alias{determineDefaultSeqnameStyle}
\alias{determineDefaultSeqnameStyle,TranscriptDb-method}
\alias{show,TranscriptDb-method}
\alias{show,FeatureDb-method}
\alias{isActiveSeq}
\alias{isActiveSeq<-}
\alias{isActiveSeq,TranscriptDb-method}
\alias{isActiveSeq<-,TranscriptDb-method}

% coercion
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
  
  See \code{?\link{saveDb}} and \code{?\link{loadDb}} for
  saving and loading the database contents of a TranscriptDb or FeatureDb
  object.
}

\section{Methods}{
  In the code snippets below, \code{x} is a TranscriptDb object. For the
  metadata and show methods, there is also support for FeatureDb objects.

  \describe{
    \item{}{
      \code{metadata(x)}:
      Return \code{x}'s metadata in a data frame.
    }
    \item{}{
      \code{seqinfo(x)}:
      Get the information about the underlying sequences as a
      \link[GenomicRanges]{Seqinfo} object.
      Note that there is no \code{seqinfo} setter for TranscriptDb
      objects (this information is stored in the SQLite db, and,
      in normal use, this db is accessed in read-only mode).
    }
    \item{}{
      \code{seqnameStyle(x)}:
      List the matching seqname styles for \code{x}.
      \code{seqnameStyle(x)} is equivalent to \code{seqnameStyle(seqinfo(x))}.
      Note that this information is not stored in \code{x} but inferred
      by looking up \code{seqlevels(x)} against a seqname style database
      stored in the seqnames.db metadata package (required).
    }
    \item{}{
      \code{determineDefaultSeqnameStyle(x)}:
      Determine the default seqname style for the database in \code{x}.
    }
    \item{}{
      \code{as.list(x)}:
      Dumps the entire db into a list of data frames \code{txdump} that can
      be used in \code{do.call(makeTranscriptDb, txdump)} to make the db
      again with no loss of information.
      Note that the transcripts are dumped in the same order in all the
      data frames.
    }
    \item{}{
      \code{isActiveSeq(x)}:
      Return the currently active sequences for this txdb object as a
      named logical vector.  Only active sequences will be tapped when
      using the supplied accessor methods.  Inactive sequences will be
      ignored.  By default, all available sequences will be active.
    }    
    \item{}{
      \code{isActiveSeq(x) <- value}:
      Allows the user to change which sequences will be actively
      accessed by the accessor methods by altering the contents of this
      named logical vector.
    }
  }

  See \code{?\link{transcripts}}, \code{?\link{transcriptsByOverlaps}},
  \code{?\link{id2name}} and \code{?\link{transcriptsBy}} for other useful
  operations on TranscriptDb objects.
}

\author{H. Pages, Marc Carlson}

\seealso{
  \itemize{
    \item \code{\link{makeTranscriptDbFromUCSC}} and
          \code{\link{makeTranscriptDbFromBiomart}} for convenient ways to
          make TranscriptDb objects from UCSC or Biomart online resources.
    \item \code{\link{transcripts}}, \code{\link{transcriptsBy}},
          and \code{\link{transcriptsByOverlaps}},
          for how to extract genomic features from a TranscriptDb object.
    \item \link[GenomicFeatures]{select-methods} for how to use the
          simple "select" interface to extract information from a
          TranscriptDb object.
    \item \code{\link{id2name}} for mapping TranscriptDb internal ids
          to external names for a given feature type.
    \item The \link[GenomicRanges]{Seqinfo} class in the GenomicRanges package.
  }
}

\examples{
txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                         package="GenomicFeatures")
txdb <- loadDb(txdb_file)
txdb

## Use of seqinfo
seqinfo(txdb)
seqlevels(txdb)  # shortcut for 'seqlevels(seqinfo(txdb))'
seqlengths(txdb)  # shortcut for 'seqlengths(seqinfo(txdb))'
isCircular(txdb)  # shortcut for 'isCircular(seqinfo(txdb))'
names(which(isCircular(txdb)))

## Examples on how to change which sequences are active
## Set chr1 and chr3 to be inactive:
isActiveSeq(txdb) <- c("1"=FALSE, "3"=FALSE)
## Set ALL of the chromsomed to be inactive
isActiveSeq(txdb)[seqlevels(txdb)] <- FALSE
## Now set only chr1 and chr5 to be active
isActiveSeq(txdb) <- c("1"=TRUE, "5"=TRUE)

## Use of as.list
txdump <- as.list(txdb)
txdump
txdb1 <- do.call(makeTranscriptDb, txdump)
stopifnot(identical(as.list(txdb1), txdump))
}

\keyword{methods}
\keyword{classes}
