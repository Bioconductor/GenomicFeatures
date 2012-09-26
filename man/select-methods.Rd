\name{select-methods}

\alias{select-methods}

\alias{cols,TranscriptDb-method}
\alias{keytypes,TranscriptDb-method}
\alias{keys,TranscriptDb-method}
\alias{select,TranscriptDb-method}

\title{Using the "select" interface on TranscriptDb objects}

\description{
  \code{select}, \code{cols} and \code{keys} can be used together to
  extract data from a \link{TranscriptDb} object.
}

\details{
  In the code snippets below, \code{x} is a \link{TranscriptDb} object.

  \describe{
    \item{}{
      \code{keytypes(x)}:
      allows the user to discover which keytypes can be passed in to
      \code{select} or \code{keys} and the \code{keytype} argument.
    }
    \item{}{
      \code{keys(x, keytype)}:
      Return keys for the database contained in the \link{TranscriptDb}
      object .  By default it will return the "TXNAME" keys for the database,
      but if used with the \code{keytype} argument, it will return the keys
      from that keytype.
    }
    \item{}{
      \code{cols(x)}:
      Show which kinds of data can be returned for the
      \link{TranscriptDb} object.
    }
    \item{}{
      \code{select(x, keys, cols, keytype)}:
      When all the appropriate arguments are specified \code{select}
      will retrieve the matching data as a data.frame based on
      parameters for selected \code{keys} and \code{cols} and
      \code{keytype} arguments.
    }
  }

}

\author{Marc Carlson}

\seealso{
  \itemize{
    \item \code{\link{transcripts}}, \code{\link{transcriptsBy}},
          and \code{\link{transcriptsByOverlaps}}, for other ways to
          extract genomic features from a \link{TranscriptDb} object.
    \item The \link{TranscriptDb} class.
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
cols(txdb)

## call select
res <- select(txdb, head(keys(txdb, "GENEID")),
              cols=c("GENEID","TXNAME"),
              keytype="GENEID")
head(res)
}

\keyword{methods}
