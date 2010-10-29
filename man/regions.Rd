\name{regions}
\Rdversion{1.1}
\alias{transcripts_deprecated}
\alias{exons_deprecated}
\alias{introns_deprecated}
\title{
  Functions that compute genomic regions of interest.
}
\description{
  Functions that compute genomic regions of interest such as promotor,
  upstream regions etc, from the genomic locations provided in a
  UCSC-style data frame.
}
\usage{
transcripts_deprecated(genes, proximal = 500, distal = 10000)
exons_deprecated(genes)
introns_deprecated(genes)
}
\arguments{
  \item{genes}{
    A UCSC-style data frame i.e. a data frame with 1 row per transcript
    and at least the following columns: \code{"name"}, \code{"chrom"},
    \code{"strand"}, \code{"txStart"}, \code{"txEnd"},
    \code{"exonCount"}, \code{"exonStarts"}, \code{"exonEnds"},
    \code{"intronStarts"} and \code{"intronEnds"}.
    A value in any of the last 4 columns must be a comma-separated list
    of integers.
    Note that unlike what UCSC does the start values here must be 1-based,
    not 0-based.
  }
  \item{proximal}{
    The number of bases on either side of TSS and 3'-end for
    the promoter and end region, respectively.
  }
  \item{distal}{
    The number of bases on either side for upstream/downstream,
    i.e. enhancer/silencer regions.
  }
}
\details{
  The assumption made for introns is that there must be more than one
  exon, and that the introns are between the end of one exon and before
  the start of the next exon.
}
\value{
  All of these functions return a
  \code{\link[IRanges:RangedData-class]{RangedData}}
  object with a \code{gene} column with the UCSC ID of the gene.
  For \code{transcripts_deprecated}, each element corresponds to a transcript, and
  there are columns for each type of region (promoter, threeprime,
  upstream, and downstream).
  For \code{exons_deprecated}, each element corresponds to an exon.
  For \code{introns_deprecated}, each element corresponds to an intron.
}
\author{
  M. Lawrence.
}

