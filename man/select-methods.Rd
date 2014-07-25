\name{select-methods}

\alias{select-methods}

\alias{columns,TxDb-method}
\alias{keytypes,TxDb-method}
\alias{keys,TxDb-method}
\alias{select,TxDb-method}

\title{Using the "select" interface on TxDb objects}

\description{
  \code{select}, \code{columns} and \code{keys} can be used together to
  extract data from a \link{TxDb} object.
}

\details{
  In the code snippets below, \code{x} is a \link{TxDb} object.

  \describe{
    \item{}{
      \code{keytypes(x)}:
      allows the user to discover which keytypes can be passed in to
      \code{select} or \code{keys} and the \code{keytype} argument.
    }
    \item{}{
      \code{keys(x, keytype, pattern, column, fuzzy)}: Return keys for
      the database contained in the \link{TxDb} object .

      The \code{keytype} argument specifies the kind of keys that will
      be returned. By default \code{keys} will return the "GENEID" keys
      for the database.

      If \code{keys} is used with \code{pattern}, it will pattern match
      on the \code{keytype}.

      But if the \code{column} argument is also provided along with the
      \code{pattern} argument, then \code{pattern} will be matched
      against the values in \code{column} instead.

      And if \code{keys} is called with \code{column} and no
      \code{pattern} argument, then it will return all keys that have
      corresponding values in the \code{column} argument.

      Thus, the behavior of \code{keys} all depends on how many arguments are
      specified.

      Use of the \code{fuzzy} argument will toggle fuzzy matching to
      TRUE or FALSE.  If \code{pattern} is not used, fuzzy is ignored.
    }
    \item{}{
      \code{columns(x)}:
      Show which kinds of data can be returned for the
      \link{TxDb} object.
    }
    \item{}{
      \code{select(x, keys, columns, keytype)}:
      When all the appropriate arguments are specified \code{select}
      will retrieve the matching data as a data.frame based on
      parameters for selected \code{keys} and \code{columns} and
      \code{keytype} arguments.
    }
  }

}



\author{Marc Carlson}

\seealso{
  \itemize{
    \item \link[AnnotationDbi]{AnnotationDb-class} for more descriptsion
    of methods \code{select},\code{keytypes},\code{keys} and \code{columns}.
    \item \code{\link{transcripts}}, \code{\link{transcriptsBy}},
          and \code{\link{transcriptsByOverlaps}}, for other ways to
          extract genomic features from a \link{TxDb} object.
    \item The \link{TxDb} class.
  }
}

\examples{
txdb_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
                         package="GenomicFeatures")
txdb <- loadDb(txdb_file)
txdb

## find key types
keytypes(txdb)

## list IDs that can be used to filter
head(keys(txdb, "GENEID"))
head(keys(txdb, "TXID"))
head(keys(txdb, "TXNAME"))

## list columns that can be returned by select
columns(txdb)

## call select
res <- select(txdb, head(keys(txdb, "GENEID")),
              columns=c("GENEID","TXNAME"),
              keytype="GENEID")
head(res)
}

\keyword{methods}
