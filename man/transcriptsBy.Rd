\name{transcriptsBy}

\alias{transcriptsBy}
\alias{transcriptsBy,TxDb-method}
\alias{exonsBy}
\alias{exonsBy,TxDb-method}
\alias{cdsBy}
\alias{cdsBy,TxDb-method}
\alias{intronsByTranscript}
\alias{intronsByTranscript,TxDb-method}
\alias{fiveUTRsByTranscript}
\alias{fiveUTRsByTranscript,TxDb-method}
\alias{threeUTRsByTranscript}
\alias{threeUTRsByTranscript,TxDb-method}

\title{
  Extract and group genomic features of a given type
}
\description{
  Generic functions to extract genomic features of a given type
  grouped based on another type of genomic feature.
  This page documents the methods for \link{TxDb} objects only.
}
\usage{
transcriptsBy(x, by=c("gene", "exon", "cds"), ...)
\S4method{transcriptsBy}{TxDb}(x, by=c("gene", "exon", "cds"), use.names=FALSE)

exonsBy(x, by=c("tx", "gene"), ...)
\S4method{exonsBy}{TxDb}(x, by=c("tx", "gene"), use.names=FALSE)

cdsBy(x, by=c("tx", "gene"), ...)
\S4method{cdsBy}{TxDb}(x, by=c("tx", "gene"), use.names=FALSE)

intronsByTranscript(x, ...)
\S4method{intronsByTranscript}{TxDb}(x, use.names=FALSE)

fiveUTRsByTranscript(x, ...)
\S4method{fiveUTRsByTranscript}{TxDb}(x, use.names=FALSE)

threeUTRsByTranscript(x, ...)
\S4method{threeUTRsByTranscript}{TxDb}(x, use.names=FALSE)
}
\arguments{
  \item{x}{A \link{TxDb} object.}
  \item{...}{Arguments to be passed to or from methods.}
  \item{by}{One of \code{"gene"}, \code{"exon"}, \code{"cds"} or \code{"tx"}.
    Determines the grouping.}
  \item{use.names}{Controls how to set the names of the returned
    \link[GenomicRanges]{GRangesList} object.
    These functions return all the features of a given type (e.g.
    all the exons) grouped by another feature type (e.g. grouped by
    transcript) in a \link[GenomicRanges]{GRangesList} object.
    By default (i.e. if \code{use.names} is \code{FALSE}), the
    names of this \link[GenomicRanges]{GRangesList} object
    (aka the group names) are the internal ids of the features
    used for grouping (aka the grouping features), which are
    guaranteed to be unique.
    If \code{use.names} is \code{TRUE}, then the names of the
    grouping features are used instead of their internal ids.
    For example, when grouping by transcript (\code{by="tx"}),
    the default group names are the transcript internal ids
    (\code{"tx_id"}). But, if \code{use.names=TRUE}, the group
    names are the transcript names (\code{"tx_name"}).
    Note that, unlike the feature ids, the feature names are not
    guaranteed to be unique or even defined (they could be all
    \code{NA}s). A warning is issued when this happens.
    See \code{?\link{id2name}} for more information about
    feature internal ids and feature external names and how
    to map the formers to the latters.

    Finally, \code{use.names=TRUE} cannot be used when grouping
    by gene \code{by="gene"}. This is because, unlike for the
    other features, the gene ids are external ids (e.g. Entrez
    Gene or Ensembl ids) so the db doesn't have a \code{"gene_name"}
    column for storing alternate gene names.
  }
}
\details{
  These functions return a \link[GenomicRanges]{GRangesList} object
  where the ranges within each of the elements are ordered according
  to the following rule:

  When using \code{exonsBy} or \code{cdsBy} with \code{by = "tx"},
  the returned exons or CDS are ordered by ascending rank for each
  transcript, that is, by their position in the transcript.
  In all other cases, the ranges will be ordered by chromosome, strand,
  start, and end values.
}
\value{A \link[GenomicRanges]{GRangesList} object.}
\author{
  M. Carlson, P. Aboyoun and H. Pagès
}
\seealso{
  \itemize{
    \item \code{\link{transcripts}} and \code{\link{transcriptsByOverlaps}}
          for more ways to extract genomic features
          from a \link{TxDb} object.

    \item \code{\link{transcriptLengths}} for extracting the transcript
          lengths from a \link{TxDb} object.

    \item \code{\link{extractTranscriptSeqs}} for extracting transcript
          (or CDS) sequences from chromosome sequences.

    \item \code{\link{coverageByTranscript}} for computing coverage by
          transcript (or CDS) of a set of ranges.

    \item \link[GenomicFeatures]{select-methods} for how to use the
          simple "select" interface to extract information from a
          \link{TxDb} object.

    \item \code{\link{id2name}} for mapping \link{TxDb} internal ids
          to external names for a given feature type.

    \item The \link{TxDb} class.
  }
}
\examples{
txdb_file <- system.file("extdata", "hg19_knownGene_sample.sqlite",
                         package="GenomicFeatures")
txdb <- loadDb(txdb_file)

## Get the transcripts grouped by gene:
transcriptsBy(txdb, "gene")

## Get the exons grouped by gene:
exonsBy(txdb, "gene")

## Get the CDS grouped by transcript:
cds_by_tx0 <- cdsBy(txdb, "tx")
## With more informative group names:
cds_by_tx1 <- cdsBy(txdb, "tx", use.names=TRUE)
## Note that 'cds_by_tx1' can also be obtained with:
names(cds_by_tx0) <- id2name(txdb, feature.type="tx")[names(cds_by_tx0)]
stopifnot(identical(cds_by_tx0, cds_by_tx1))

## Get the introns grouped by transcript:
intronsByTranscript(txdb)

## Get the 5' UTRs grouped by transcript:
fiveUTRsByTranscript(txdb)
fiveUTRsByTranscript(txdb, use.names=TRUE)  # more informative group names
}

\keyword{methods}
