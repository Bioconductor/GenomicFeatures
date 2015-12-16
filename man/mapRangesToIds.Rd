\docType{methods}
\name{mapRangesToIds}
\alias{mapRangesToIds}
\alias{mapRangesToIds,TxDb-method}
\title{Map Genomic Ranges to IDs}
\usage{
mapRangesToIds(x, ...)

\S4method{mapRangesToIds}{TxDb}(x, ranges, type = c("cds", "exon", "tx",
  "gene"), columns = NULL, ...)
}
\arguments{
\item{x}{Database to use for mapping}

\item{ranges}{range object used to subset}

\item{type}{of feature to return}

\item{columns}{additional metadata columns to include in the output.}

\item{...}{Additional arguments passed to
  \code{\link[GenomicRanges:findOverlaps-methods]{findOverlaps}}}
}
\value{
\code{\link[S4Vectors]{DataFrame}} of mcols from the database.
}
\description{
Map Genomic Ranges to IDs
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
res <- mapIdsToRanges(txdb, keys = keys, type = "tx")
mapRangesToIds(txdb, res, "tx")
}
