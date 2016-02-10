\name{TxDb-class}

\alias{TxDb-class}
\alias{class:TxDb}
\alias{TxDb}

\alias{species,TxDb-method}
\alias{organism,TxDb-method}
\alias{seqlevels0,TxDb-method}
\alias{seqlevels<-,TxDb-method}
\alias{seqinfo,TxDb-method}
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
      \code{seqlevels0(x)}:
      Get the \emph{sequence levels} originally in \code{x}. This ignores any
      change the user might have made to the \emph{sequence levels} with the
      \code{seqlevels} setter.
    }
    \item{}{
      \code{seqlevels(x)}, \code{seqlevels(x) <- value}:
      Get or set the \emph{sequence levels} in \code{x}.
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

\author{Hervé Pagès, Marc Carlson}

\seealso{
  \itemize{
    \item \code{\link{makeTxDbFromUCSC}}, \code{\link{makeTxDbFromBiomart}},
          \code{\link{makeTxDbFromGRanges}}, and \code{\link{makeTxDbFromGFF}},
          for convenient ways to make a \link{TxDb} object from UCSC or BioMart
          online resources, or from a \link[GenomicRanges]{GRanges} object,
          or from a GFF or GTF file.

    \item \code{\link[AnnotationDbi]{saveDb}} and
          \code{\link[AnnotationDbi]{loadDb}} in the \pkg{AnnotationDbi}
          package for saving and loading a TxDb object as an SQLite file.

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

## You can set user-supplied seqlevels on 'txdb' to restrict any further
## operations to a subset of chromosomes:
seqlevels(txdb) <- c("Y", "6")
## Then you can restore the seqlevels stored in the db:
seqlevels(txdb) <- seqlevels0(txdb)

## Use of as.list():
txdb_dump <- as.list(txdb)
txdb_dump
txdb1 <- do.call(makeTxDb, txdb_dump)
stopifnot(identical(as.list(txdb1), txdb_dump))
}

\keyword{methods}
\keyword{classes}
