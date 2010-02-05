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
  upstream regions etc, from the genomic locations provided in data like
  the data.frame returned by \code{\link[GenomicFeatures.Hsapiens.UCSC.hg18]{geneHuman}}.
}
\usage{
transcripts_deprecated(genes, proximal = 500, distal = 10000)
exons_deprecated(genes)
introns_deprecated(genes)
}
\arguments{
  \item{genes}{ A data.frame like that returned by
    \code{\link[GenomicFeatures.Hsapiens.UCSC.hg18]{geneHuman}}.}
  \item{proximal}{
    The number of bases on either side of TSS and 3'-end for
    the promoter and end region, respectively.}
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

\examples{
library(GenomicFeatures.Hsapiens.UCSC.hg18)
## promoter 300bp up and down from TSS (threeprime from TES)
transcripts_deprecated(geneHuman(), proximal = 300)
}
