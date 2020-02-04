\name{makeTxDbFromEnsembl}

\alias{makeTxDbFromEnsembl}

\title{
  Make a TxDb object from an Ensembl database
}

\description{
  The \code{makeTxDbFromEnsembl} function creates a \link{TxDb} object for
  a given organism by importing the genomic locations of its transcripts,
  exons, CDS, and genes from an Ensembl database.

  Note that it uses the \pkg{RMariaDB} package internally so make sure
  that this package is installed.
}

\usage{
makeTxDbFromEnsembl(organism="Homo sapiens",
                    release=NA,
                    circ_seqs=NULL,
                    server="ensembldb.ensembl.org",
                    username="anonymous", password=NULL, port=0L,
                    tx_attrib=NULL)
}

\arguments{
  \item{organism}{
    The \emph{scientific name} (i.e. genus and species, or genus and species
    and subspecies) of the organism for which to import the data.
    Case is not sensitive. Underscores can be used instead of white spaces
    e.g. \code{"homo_sapiens"} is accepted.
  }
  \item{release}{
    The Ensembl release to query e.g. 89. If set to \code{NA} (the default),
    the current release is used.
  }
  \item{circ_seqs}{
    A character vector to list out which chromosomes should be marked
    as circular.
  }
  \item{server}{
    The name of the MySQL server to query.
    See \url{https://www.ensembl.org/info/data/mysql.html} for the list of
    Ensembl public MySQL servers.
    Make sure to use the server nearest to you. It can make a big difference!
  }
  \item{username}{
    Login username for the MySQL server.
  }
  \item{password}{
    Login password for the MySQL server.
  }
  \item{port}{
    Port of the MySQL server.
  }
  \item{tx_attrib}{
    If not \code{NULL}, only select transcripts with an attribute of
    the given code, a string, like \code{"gencode_basic"}.
  }
}

\value{
  A \link{TxDb} object.
}

\note{
  \code{makeTxDbFromEnsembl} tends to be faster and more reliable than
  \code{\link{makeTxDbFromBiomart}}.
}

\author{H. Pagès}

\seealso{
  \itemize{
    \item \code{\link{makeTxDbFromUCSC}} and \code{\link{makeTxDbFromBiomart}}
          for making a \link{TxDb} object from other online resources.

    \item \code{\link{makeTxDbFromGRanges}} and \code{\link{makeTxDbFromGFF}}
          for making a \link{TxDb} object from a \link[GenomicRanges]{GRanges}
          object, or from a GFF or GTF file.

    \item The \link{TxDb} class.

    \item \code{\link{makeTxDb}} for the low-level function used by the
          \code{makeTxDbFrom*} functions to make the \link{TxDb} object
          returned to the user.
  }
}

\examples{
\dontrun{
txdb <- makeTxDbFromEnsembl("Saccharomyces cerevisiae",
                            server="useastdb.ensembl.org")
txdb
}
}
