\docType{methods}
\name{mapIdsToRanges}
\alias{mapIdsToRanges}
\alias{mapIdsToRanges,TxDb-method}
\title{Map IDs to Genomic Ranges}
\usage{
mapIdsToRanges(db, ...)

\S4method{mapIdsToRanges}{TxDb}(db, vals, type = c("cds", "exon", "tx",
  "gene"), columns = NULL, ...)
}
\arguments{
\item{db}{Database to use for mapping}

\item{...}{Additional arguments passed to methods}

\item{vals}{Values to lookup, passed to \code{\link{transcripts}} et. al.}

\item{type}{Types of feature to return}

\item{columns}{Additional metadata columns to include in the output}
}
\value{
\code{\link[GenomicRanges]{GRangesList}} corresponding to the vals
}
\description{
Map IDs to Genomic Ranges
}
\section{Methods (by class)}{
\itemize{
\item \code{TxDb}: TxDb method
}}
\examples{
txdb <- makeTxDbFromGRanges(
    readRDS(system.file(package = "GenomicFeatures", "extdata", "sample_ranges.rds")))

vals <- list(tx_name = c("ENST00000371582", "ENST00000371588", "ENST00000494752", "ENST00000614008", "ENST00000496771"))
mapIdsToRanges(txdb, vals = vals, type = "tx")
}
