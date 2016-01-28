\name{transcripts}

\alias{transcripts}
\alias{transcripts,TxDb-method}
\alias{exons}
\alias{exons,TxDb-method}
\alias{cds}
\alias{cds,TxDb-method}
\alias{genes}
\alias{genes,TxDb-method}
\alias{promoters}
\alias{promoters,TxDb-method}

\title{
  Extract genomic features from an object
}

\description{
  Generic functions to extract genomic features from an object.
  This page documents the methods for \link{TxDb} objects only.
}

\usage{
transcripts(x, ...)
\S4method{transcripts}{TxDb}(x, vals=NULL, columns=c("tx_id", "tx_name"))

exons(x, ...)
\S4method{exons}{TxDb}(x, vals=NULL, columns="exon_id")

cds(x, ...)
\S4method{cds}{TxDb}(x, vals=NULL, columns="cds_id")

genes(x, ...)
\S4method{genes}{TxDb}(x, vals=NULL, columns="gene_id", single.strand.genes.only=TRUE)

\S4method{promoters}{TxDb}(x, upstream=2000, downstream=200, ...)
}

\arguments{ 
  \item{x}{
    A \link{TxDb} object.
  }
  \item{...}{
    For the \code{transcripts}, \code{exons}, \code{cds}, and \code{genes}
    generic functions: arguments to be passed to methods.

    For the \code{promoters} method for \link{TxDb} objects: arguments
    to be passed to the internal call to \code{transcripts}.
  }
  \item{vals}{
    Either \code{NULL} or a named list of vectors to be used to
    restrict the output. Valid names for this list are:
    \code{"gene_id"},
    \code{"tx_id"}, \code{"tx_name"},
    \code{"tx_chrom"}, \code{"tx_strand"},
    \code{"exon_id"}, \code{"exon_name"},
    \code{"exon_chrom"}, \code{"exon_strand"},
    \code{"cds_id"}, \code{"cds_name"},
    \code{"cds_chrom"}, \code{"cds_strand"} and \code{"exon_rank"}.
  }
  \item{columns}{
    Columns to include in the output.
    Must be \code{NULL} or a character vector as given by the
    \code{columns} method. With the following restrictions:
    \itemize{
      \item \code{"TXCHROM"} and \code{"TXSTRAND"} are not allowed
      for \code{transcripts}.
      \item \code{"EXONCHROM"} and \code{"EXONSTRAND"} are not allowed
      for \code{exons}.
      \item \code{"CDSCHROM"} and \code{"CDSSTRAND"} are not allowed
      for \code{cds}.
    }
    If the vector is named, those names are used for the corresponding
    column in the element metadata of the returned object. 
  }
  \item{single.strand.genes.only}{
    \code{TRUE} or \code{FALSE}. If \code{TRUE} (the default), then genes
    that have exons located on both strands of the same chromosome or on
    two different chromosomes are dropped. In that case, the genes are
    returned in a \link[GenomicRanges]{GRanges} object. Otherwise, all genes
    are returned in a \link[GenomicRanges]{GRangesList} object with the
    columns specified thru the \code{columns} argument set as \emph{top level}
    metadata columns. (Please keep in mind that the \emph{top level} metadata
    columns of a \link[GenomicRanges]{GRangesList} object are not displayed
    by the \code{show} method.)
  }
  \item{upstream}{
    For \code{promoters} : An \code{integer(1)} value indicating the
    number of bases upstream from the transcription start site. For
    additional details see \code{?`promoters,GRanges-method`}.
  }
  \item{downstream}{
    For \code{promoters} : An \code{integer(1)} value indicating the 
    number of bases downstream from the transcription start site. For 
    additional details see \code{?`promoters,GRanges-method`}.
  }
}

\details{
  These are the main functions for extracting transcript information
  from a \link{TxDb} object. These methods can restrict the output based
  on categorical information. To restrict the output based on interval
  information, use the \code{\link{transcriptsByOverlaps}},
  \code{\link{exonsByOverlaps}}, and \code{\link{cdsByOverlaps}}
  functions.

  The \code{promoters} function computes user-defined promoter regions
  for the transcripts in a \link{TxDb} object. The return object is a 
  \code{GRanges} of promoter regions around the transcription start 
  site the span of which is defined by \code{upstream} and \code{downstream}.
  For additional details on how the promoter range is computed and the
  handling of \code{+} and \code{-} strands see 
  \code{?`promoters,GRanges-method`}.
}

\value{
  A \link[GenomicRanges]{GRanges} object. The only exception being
  when \code{genes} is used with \code{single.strand.genes.only=FALSE},
  in which case a \link[GenomicRanges]{GRangesList} object is returned.
}

\author{
  M. Carlson, P. Aboyoun and H. Pages
}

\seealso{
  \itemize{
    \item \code{\link{transcriptsBy}} and \code{\link{transcriptsByOverlaps}}
          for more ways to extract genomic features
          from a \link{TxDb} object.

    \item \code{\link{transcriptLengths}} for extracting the transcript
          lengths from a \link{TxDb} object.

    \item \code{\link{extractTranscriptSeqs}} for extracting transcript
          (or CDS) sequences from chromosome sequences.

    \item \code{\link{coverageByTranscript}} for computing coverage by
          transcript (or CDS) of a set of ranges.

    \item \code{\link{disjointExons}} for extracting the non-overlapping
          exon parts from a \link{TxDb} object.

    \item \link[GenomicFeatures]{select-methods} for how to use the
          simple "select" interface to extract information from a
          \link{TxDb} object.

    \item \code{\link{microRNAs}} and \code{\link{tRNAs}} for extracting
          microRNA or tRNA genomic ranges from a \link{TxDb} object.

    \item \code{\link{id2name}} for mapping \link{TxDb} internal ids
          to external names for a given feature type.

    \item The \link{TxDb} class.
  }
}

\examples{
txdb_file <- system.file("extdata", "hg19_knownGene_sample.sqlite",
                         package="GenomicFeatures")
txdb <- loadDb(txdb_file)

## ---------------------------------------------------------------------
## transcripts()
## ---------------------------------------------------------------------

tx <- transcripts(txdb)
tx

## A sanity check:
stopifnot(identical(mcols(tx)$tx_id, seq_along(tx)))

vals <- list(tx_chrom = c("chr3", "chr5"), tx_strand = "+")
transcripts(txdb, vals)

## ---------------------------------------------------------------------
## exons()
## ---------------------------------------------------------------------

exons(txdb, vals=list(exon_id=1), columns=c("EXONID", "TXNAME"))
exons(txdb, vals=list(tx_name="uc009vip.1"), columns=c("EXONID",
      "TXNAME"))

## ---------------------------------------------------------------------
## genes()
## ---------------------------------------------------------------------

genes(txdb)  # a GRanges object
cols <- c("tx_id", "tx_chrom", "tx_strand",
          "exon_id", "exon_chrom", "exon_strand")
single_strand_genes <- genes(txdb, columns=cols)

## Because we've returned single strand genes only, the "tx_chrom"
## and "exon_chrom" metadata columns are guaranteed to match
## 'seqnames(single_strand_genes)':
stopifnot(identical(as.character(seqnames(single_strand_genes)),
                    as.character(mcols(single_strand_genes)$tx_chrom)))
stopifnot(identical(as.character(seqnames(single_strand_genes)),
                    as.character(mcols(single_strand_genes)$exon_chrom)))
## and also the "tx_strand" and "exon_strand" metadata columns are
## guaranteed to match 'strand(single_strand_genes)':
stopifnot(identical(as.character(strand(single_strand_genes)),
                    as.character(mcols(single_strand_genes)$tx_strand)))
stopifnot(identical(as.character(strand(single_strand_genes)),
                    as.character(mcols(single_strand_genes)$exon_strand)))

all_genes <- genes(txdb, columns=cols, single.strand.genes.only=FALSE)
all_genes  # a GRangesList object
multiple_strand_genes <- all_genes[elementNROWS(all_genes) >= 2]
multiple_strand_genes
mcols(multiple_strand_genes)

## ---------------------------------------------------------------------
## promoters()
## ---------------------------------------------------------------------

## This:
promoters(txdb, upstream=100, downstream=50)
## is equivalent to:
promoters(transcripts(txdb), upstream=100, downstream=50)

## Extra arguments are passed to transcripts(). So this:
promoters(txdb, upstream=100, downstream=50,
          columns=c("tx_name", "gene_id"))
## is equivalent to:
promoters(transcripts(txdb, columns=c("tx_name", "gene_id")),
          upstream=100, downstream=50)
}

\keyword{methods}
