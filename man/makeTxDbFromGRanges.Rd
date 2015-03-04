\name{makeTxDbFromGRanges}

\alias{makeTxDbFromGRanges}

\title{Make a TxDb object from a GRanges object}

\description{
  The \code{makeTxDbFromGRanges} function allows the user
  to extract gene, transcript, exon, and CDS information from a
  \link[GenomicRanges]{GRanges} object structured as GFF3 or GTF,
  and to return that information in a \link{TxDb} object.
}

\usage{
makeTxDbFromGRanges(gr, drop.stop.codons=FALSE, metadata=NULL)
}

\arguments{
  \item{gr}{
    A \link[GenomicRanges]{GRanges} object structured as GFF3 or GTF,
    typically obtained with \code{rtracklayer::\link[rtracklayer]{import}()}.
  }
  \item{drop.stop.codons}{
    \code{TRUE} or \code{FALSE}. If \code{TRUE}, then features of type
    \code{stop_codon} are ignored. Otherwise (the default) the stop codons
    are considered to be part of the CDS and merged to them.
  }
  \item{metadata}{
    A 2-column data frame containing meta information to be included in the
    \link{TxDb} object. This data frame is just passed to
    \code{\link{makeTxDb}}, which \code{makeTxDbFromGRanges} calls at the
    end to make the \link{TxDb} object from the information extracted from
    \code{gr}. See \code{?\link{makeTxDb}} for more information about the
    format of \code{metadata}.
  }
}

\details{
}

\value{A \link{TxDb} object.}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item \code{\link{makeTxDbFromUCSC}}, \code{\link{makeTxDbFromBiomart}},
          and \code{\link{makeTxDbFromGFF}}, for convenient ways to make a
          \link{TxDb} object from UCSC or BioMart online resources, or
          directly from a GFF or GTF file.

    \item The \code{asGFF} method for \link{TxDb} objects
          (\link{asGFF,TxDb-method}) for the reverse of
          \code{makeTxDbFromGRanges}, that is, for turning a \link{TxDb}
          object into a \link[GenomicRanges]{GRanges} object structured as
          GFF.

    \item \code{\link{makeTxDb}} for the low-level function used by the
          above functions for making the \link{TxDb} object returned to the
          user.
  }
}

\examples{
## ---------------------------------------------------------------------
## WITH A GRanges OBJECT STRUCTURED AS GFF3
## ---------------------------------------------------------------------
GFF3_files <- system.file("extdata", "GFF3_files",
                          package="GenomicFeatures")

path <- file.path(GFF3_files, "a.gff3")
gr <- import(path)
txdb <- makeTxDbFromGRanges(gr)
txdb

## Reverse operation:
gr2 <- asGFF(txdb)

## Sanity check:
stopifnot(identical(as.list(txdb), as.list(makeTxDbFromGRanges(gr2))))

## ---------------------------------------------------------------------
## WITH A GRanges OBJECT STRUCTURED AS GTF
## ---------------------------------------------------------------------

}
