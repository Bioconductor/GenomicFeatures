\name{TxDb-class}

\alias{TxDb-class}
\alias{class:TxDb}
\alias{TxDb}

\alias{seqinfo,TxDb-method}
\alias{seqinfo<-,TxDb-method}
\alias{isActiveSeq}
\alias{isActiveSeq<-}
\alias{isActiveSeq,TxDb-method}
\alias{isActiveSeq<-,TxDb-method}
\alias{show,TxDb-method}

% coercion
\alias{as.list,TxDb-method}

\title{TxDb objects}

\description{
  The TxDb class is a container for storing transcript annotations.

  See \code{?\link{FeatureDb}} for a more generic container for storing
  genomic locations of an arbitrary type of genomic features.

  See \code{?\link{makeTxDbFromUCSC}} and
  \code{?\link{makeTxDbFromBiomart}} for convenient ways to
  make TxDb objects from UCSC or BioMart online resources.

  See \code{?\link{makeTxDbFromGFF}} for making a TxDb
  object from annotations available as a GFF3 or GTF file.
}

\section{Methods}{
  In the code snippets below, \code{x} is a TxDb object.

  \describe{
    \item{}{
      \code{metadata(x)}:
      Return \code{x}'s metadata in a data frame.
    }
    \item{}{
      \code{seqinfo(x)}, \code{seqinfo(x) <- value}:
      Get or set the information about the underlying sequences.
      Note that, for now, the setter only supports replacement of the
      sequence names, i.e., except for their sequence names (accessed with
      \code{seqnames(value)} and \code{seqnames(seqinfo(x))}, respectively),
      \link[GenomeInfoDb]{Seqinfo} objects \code{value} (supplied) and
      \code{seqinfo(x)} (current) must be identical.
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
    \item{}{
      \code{seqlevelsStyle(x)}, \code{seqlevelsStyle(x) <- value}:
      Get or set the seqname style for \code{x}.
      See the \link[GenomeInfoDb]{seqlevelsStyle} generic getter and setter
      in the \pkg{GenomeInfoDb} package for more information.
    }
    \item{}{
      \code{as.list(x)}:
      Dump the entire db into a list of data frames, say \code{txdb_dump},
      that can then be used to recreate the original db with
      \code{do.call(makeTxDb, txdb_dump)} with no loss of information
      (except possibly for some of the metadata).
      Note that the transcripts are dumped in the same order in all the
      data frames.
    }
  }
}

\author{H. Pages, Marc Carlson}

\seealso{
  \itemize{
    \item \code{\link{makeTxDbFromUCSC}}, \code{\link{makeTxDbFromBiomart}},
          \code{\link{makeTxDbFromGRanges}}, and \code{\link{makeTxDbFromGFF}},
          for convenient ways to make a \link{TxDb} object from UCSC or BioMart
          online resources, or from a \link[GenomicRanges]{GRanges} object,
          or from a GFF or GTF file.

    \item \code{\link{saveDb}} and \code{\link{loadDb}} for
          saving and loading the database content of a TxDb object.

    \item \code{\link{transcripts}}, \code{\link{transcriptsBy}},
          and \code{\link{transcriptsByOverlaps}},
          for how to extract genomic features from a TxDb object.

    \item \code{\link{transcriptLengths}} for extracting the transcript
          lengths from a \link{TxDb} object.

    \item \link[GenomicFeatures]{select-methods} for how to use the
          simple "select" interface to extract information from a
          TxDb object.

    \item The \link{FeatureDb} class for storing genomic locations
          of an arbitrary type of genomic features.

    \item The \link[GenomeInfoDb]{Seqinfo} class in the \pkg{GenomeInfoDb}
          package.
  }
}

\examples{
txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                         package="GenomicFeatures")
txdb <- loadDb(txdb_file)
txdb

## Use of seqinfo():
seqlevelsStyle(txdb)
seqinfo(txdb)
seqlevels(txdb)
seqlengths(txdb)  # shortcut for 'seqlengths(seqinfo(txdb))'
isCircular(txdb)  # shortcut for 'isCircular(seqinfo(txdb))'
names(which(isCircular(txdb)))

## Can set txdb so that only chr1 and chr2 are used by using the
## seqlevels:
seqlevels(txdb, force=TRUE) <- c("1", "2")
## Then you can restore the default seqlevels:
txdb <- restoreSeqlevels(txdb)


## Use of as.list():
txdb_dump <- as.list(txdb)
txdb_dump
txdb1 <- do.call(makeTxDb, txdb_dump)
stopifnot(identical(as.list(txdb1), txdb_dump))
}

\keyword{methods}
\keyword{classes}
