\docType{methods}
\name{mapRangesToIds}
\alias{mapRangesToIds}
\alias{mapRangesToIds,TxDb-method}
\title{Map Genomic Ranges to IDs}
\usage{
mapRangesToIds(db, ...)

\S4method{mapRangesToIds}{TxDb}(db, ranges, type = c("cds", "exon", "tx",
  "gene"), columns = NULL, ...)
}
\arguments{
\item{db}{Database to use for mapping}

\item{...}{Additional arguments passed to methods}

\item{ranges}{range object used to subset}

\item{type}{of feature to return}

\item{columns}{additional metadata columns to include in the output.}
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
txdb <- makeTxDbFromGRanges(
    readRDS(system.file(package = "GenomicFeatures", "extdata", "sample_ranges.rds")))

vals <- list(tx_name = c("ENST00000371582", "ENST00000371588", "ENST00000494752", "ENST00000614008", "ENST00000496771"))
res <- mapIdsToRanges(txdb, vals = vals, type = "tx")
mapRangesToIds(txdb, res, "tx")
}
