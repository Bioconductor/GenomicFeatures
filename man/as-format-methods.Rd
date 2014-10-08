\name{as-format-methods}
\alias{asBED,TxDb-method}
\alias{asGFF,TxDb-method}

\title{Coerce to file format structures}
\description{
  These functions coerce a \code{\linkS4class{TxDb}} object to a
  \code{\link[GenomicRanges:GRanges-class]{GRanges}} object with
  metadata columns encoding transcript structures according to the
  model of a standard file format. Currently, BED and GFF models are
  supported. If a \code{TxDb} is passed to
  \code{\link[rtracklayer]{export}}, when targeting a BED or GFF file,
  this coercion occurs automatically.
}
\usage{
\S4method{asBED}{TxDb}(x)
\S4method{asGFF}{TxDb}(x)
}

\arguments{
  \item{x}{
    A \code{TxDb} object to coerce to a \code{GRanges},
    structured as BED or GFF.
  }
}

\value{
  For \code{asBED}, a \code{GRanges}, with the columns \code{name},
  \code{thickStart}, \code{thickEnd}, \code{blockStarts},
  \code{blockSizes} added. The thick regions correspond to the CDS
  regions, and the blocks represent the exons. The transcript IDs are
  stored in the \code{name} column. The ranges are the transcript bounds.

  For \code{asGFF}, a \code{GRanges}, with columns \code{type},
  \code{Name}, \code{ID},, and \code{Parent}. The gene structures are
  expressed according to the conventions defined by the GFF3 spec. There
  are elements of each \code{type} of feature: \dQuote{gene},
  \dQuote{mRNA} \dQuote{exon} and \dQuote{cds}. The \code{Name} column
  contains the \code{gene_id} for genes, \code{tx_name} for transcripts,
  and exons and cds regions are \code{NA}. The \code{ID} column uses
  \code{gene_id} and \code{tx_id}, with the prefixes \dQuote{GeneID} and
  \dQuote{TxID} to ensure uniqueness across types. The exons and cds
  regions have \code{NA} for \code{ID}. The \code{Parent} column
  contains the \code{ID}s of the parent features. A feature may have
  multiple parents (the column is a \code{CharacterList}). Each exon
  belongs to one or more mRNAs, and mRNAs belong to a gene.
}

\author{
  Michael Lawrence
}

\examples{
  txdb_file <- system.file("extdata", "hg19_knownGene_sample.sqlite",
                           package="GenomicFeatures")
  txdb <- loadDb(txdb_file)

  asBED(txdb)
  asGFF(txdb)
}
