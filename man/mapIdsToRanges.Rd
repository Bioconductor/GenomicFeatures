\docType{methods}
\name{mapIdsToRanges}
\alias{mapIdsToRanges}
\alias{mapIdsToRanges,TxDb-method}
\title{Map IDs to Genomic Ranges}
\usage{
mapIdsToRanges(x, ...)

\S4method{mapIdsToRanges}{TxDb}(x, keys, type = c("cds", "exon", "tx",
  "gene"), columns = NULL)
}
\arguments{
\item{x}{Database to use for mapping}

\item{keys}{Values to lookup, passed to \code{\link{transcripts}} et. al.}

\item{type}{Types of feature to return}

\item{columns}{Additional metadata columns to include in the output}

\item{...}{Additional arguments passed to methods}
}
\value{
\code{\link[GenomicRanges]{GRangesList}} corresponding to the keys
}
\description{
Map IDs to Genomic Ranges
}
\section{Methods (by class)}{
\itemize{
\item \code{TxDb}: TxDb method
}}
\examples{
fl <- system.file(package = "GenomicFeatures", "extdata", "sample_ranges.rds")
txdb <- makeTxDbFromGRanges(readRDS(fl))

keys <- list(tx_name = c("ENST00000371582", "ENST00000371588",
    "ENST00000494752", "ENST00000614008", "ENST00000496771"))
mapIdsToRanges(txdb, keys = keys, type = "tx")
}
