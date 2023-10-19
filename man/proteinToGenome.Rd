\name{proteinToGenome}

\alias{proteinToGenome}
\alias{proteinToGenome,GRangesList-method}
\alias{proteinToGenome,ANY-method}

\title{Map protein-relative coordinates to genomic coordinates}

\description{
  \code{proteinToGenome} is a generic function for mapping
  ranges of protein-relative positions to the genome.

  NOTE: This man page is for the \code{proteinToGenome} S4 generic
  function and methods defined in the \pkg{GenomicFeatures} package,
  which are (loosely) modeled on the \code{\link[ensembldb]{proteinToGenome}}
  function from the \pkg{ensembldb} package.
  See \code{?ensembldb::\link[ensembldb]{proteinToGenome}} for the latter.
}

\usage{
## S4 generic function:
proteinToGenome(x, db, ...)  # dispatch is on 2nd argument 'db'

\S4method{proteinToGenome}{ANY}(x, db)

\S4method{proteinToGenome}{GRangesList}(x, db)
}

\arguments{
  \item{x}{
    A named \link[IRanges]{IRanges} object (or derivative) containing ranges
    of \emph{protein-relative positions} (protein-relative positions are
    positions relative to a protein sequence).

    The names on \code{x} must be transcript names present in \code{db}.
    More precisely, for the default \code{proteinToGenome()} method,
    \code{names(x)} must be a subset of:
    \preformatted{  mcols(transcripts(db, columns="tx_name"))$tx_name
    }
    And for the method for \link[GenomicRanges]{GRangesList} objects,
    \code{names(x)} must be a subset of:
    \preformatted{  names(db)
    }
  }
  \item{db}{
    For the default \code{proteinToGenome()} method: A \link{TxDb}
    object or any object that supports \code{\link{transcripts}()}
    and \code{\link{cdsBy}()} (e.g. an \link[ensembldb]{EnsDb} object
    from the \pkg{ensembldb} package).

    For the method for \link[GenomicRanges]{GRangesList} objects:
    A named \link[GenomicRanges]{GRangesList} object (or derivative)
    where each list element is a \link[GenomicRanges]{GRanges} object
    representing a CDS (the ranges in the \link[GenomicRanges]{GRanges}
    object must represent the CDS parts ordered by ascending exon rank).
  }
  \item{...}{
    Further arguments to be passed to specific methods.
  }
}

\details{
  The \code{proteinToGenome()} method for \link[GenomicRanges]{GRangesList}
  objects is the workhorse behind the default method. Note that the latter
  is a thin wrapper around the former, which simply does the following:
  \enumerate{
    \item Use \code{\link{cdsBy}()} to extract the CDS parts from \code{db}.
          The CDS parts are returned in a \link[GenomicRanges]{GRangesList}
          object that has the names of the transcript on it (one transcript
          name per list element).
    \item Call \code{proteinToGenome()} on \code{x} and the
          \link[GenomicRanges]{GRangesList} object returned by
          \code{\link{cdsBy}()}.
  }
}

\value{
  A named \link[GenomicRanges]{GRangesList} object \emph{parallel} to
  \code{x} (the transcript names on \code{x} are propagated).
  The i-th list element in the returned object is the result of mapping
  the range of protein-relative positions \code{x[i]} to the genome.

  Note that a given range in \code{x} can only be mapped to the genome
  if the name on it is the name of a \emph{coding} transcript. If it's
  not (i.e. if it's the name of a \emph{non-coding} transcript), then
  an empty \link[GenomicRanges]{GRanges} object is placed in the returned
  object to indicate the impossible mapping, and a warning is issued.

  Otherwise, if a given range in \code{x} can be mapped to the
  genome, then the result of the mapping is represented by a
  non-empty \link[GenomicRanges]{GRanges} object.
  Note that this object represents the original CDS associated to
  \code{x}, trimmed on its 5' end or 3' end, or on both.
  Furthermore, this object will have the same metadata columns as the
  \link[GenomicRanges]{GRanges} object representing the original CDS,
  plus the 2 following ones:
  \itemize{
    \item \code{protein_start}: The protein-relative start of the mapping.
    \item \code{protein_end}: The protein-relative end of the mapping.
  }
}

\note{
  Unlike \code{ensembldb::\link[ensembldb]{proteinToGenome}()} which
  can work either with Ensembl protein IDs or Ensembl transcript IDs
  on \code{x}, the default \code{proteinToGenome()} method described
  above only accepts \emph{transcript names} on \code{x}.

  This means that, if the user is in possession of protein IDs, they
  must first replace them with the corresponding transcript IDs (referred
  to as \emph{transcript names} in the context of \link{TxDb} objects).
  How to do this exactly depends on the origin of those IDs (UCSC,
  Ensembl, GTF/GFF3 file, FlyBase, etc...)
}

\author{H. Pagès, using \code{ensembldb::proteinToGenome()} for
        inspiration and design.}

\seealso{
  \itemize{
    \item The \code{\link[ensembldb]{proteinToGenome}} function in the
          \pkg{ensembldb} package, which the \code{proteinToGenome()}
          generic and methods documented in this man page are (loosely)
          modeled on.

    \item \link{TxDb} objects.

    \item \link[ensembldb]{EnsDb} objects (\link{TxDb}-like objects) in
          the \pkg{ensembldb} package.

    \item \code{\link{transcripts}} for extracting transcripts from a
          \link{TxDb}-like object.

    \item \code{\link{cdsBy}} for extracting CDS parts from a
          \link{TxDb}-like object.

    \item \link[IRanges]{IRanges} objects in the \pkg{IRanges} package.

    \item \link[GenomicRanges]{GRanges} and \link[GenomicRanges]{GRangesList}
          objects in the \pkg{GenomicRanges} package.
  }
}

\examples{
## ---------------------------------------------------------------------
## USING TOY CDS
## ---------------------------------------------------------------------

## CDS1 has 2 CDS parts:
CDS1 <- GRanges(c("chrX:11-60:+", "chrX:101-125:+"))

## CDS2 has 3 CDS parts:
CDS2 <- GRanges(c("chrY:201-230:-", "chrY:101-125:-", "chrY:11-60:-"))

## Put them in a GRangesList object:
cds_by_tx <- GRangesList(TX1=CDS1, TX2=CDS2)
cds_by_tx

x1 <- IRanges(start=8, end=20, names="TX1")
proteinToGenome(x1, cds_by_tx)

x2 <- IRanges(start=c(1, 18), end=c(25, 20), names=c("TX1", "TX1"))
x2
proteinToGenome(x2, cds_by_tx)

x3 <- IRanges(start=8, end=15, names="TX2")
proteinToGenome(x3, cds_by_tx)

x4 <- c(x3, x2)
x4
proteinToGenome(x4, cds_by_tx)

## ---------------------------------------------------------------------
## USING A TxDb OBJECT
## ---------------------------------------------------------------------
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene

## The first transcript (FBtr0309810) is non-coding:
x <- IRanges(c(FBtr0309810="11-55", FBtr0306539="90-300"))
res <- proteinToGenome(x, txdb)
res
}

\keyword{methods}
\keyword{utilities}
