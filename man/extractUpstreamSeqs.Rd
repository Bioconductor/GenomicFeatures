\name{extractUpstreamSeqs}

\alias{extractUpstreamSeqs}
\alias{extractUpstreamSeqs,GenomicRanges-method}
\alias{extractUpstreamSeqs,TranscriptDb-method}
\alias{extractUpstreamSeqs,GRangesList-method}


\title{Extract sequences upstream of a set of genes or transcripts}

\description{
  \code{extractUpstreamSeqs} is a generic function for extracting
  sequences upstream of a supplied set of genes or transcripts.
}

\usage{
extractUpstreamSeqs(x, genes, width=1000, ...)

## Dispatch is on the 2nd argument!

\S4method{extractUpstreamSeqs}{GenomicRanges}(x, genes, width=1000)

\S4method{extractUpstreamSeqs}{TranscriptDb}(x, genes, width=1000, exclude.seqlevels=NULL)
}

\arguments{
  \item{x}{
    A \link[BSgenome]{BSgenome}, \link[rtracklayer]{TwoBitFile}, or
    \link[Rsamtools]{FaFile} object, or any object containing the chromosome
    sequences from which to extract the upstream sequences. More formally,
    \code{x} must be an object for which \code{\link[GenomeInfoDb]{seqinfo}}
    and \code{\link[Biostrings]{getSeq}} are defined.
  }
  \item{genes}{
    An object containing the locations (i.e. chromosome name, start, end, and
    strand) of the genes or transcripts with respect to the reference genome.
    Only \link[GenomicRanges]{GenomicRanges} and \link{TranscriptDb} objects
    are supported at the moment. If the latter, the gene locations are obtained
    by calling \code{\link{genes}} on the \link{TranscriptDb} object
    internally.
  }
  \item{width}{
    How many bases to extract upstream of each TSS (transcription start site).
  }
  \item{...}{
    Additional arguments, for use in specific methods.
  }
  \item{exclude.seqlevels}{
    A character vector containing the chromosome names (a.k.a. sequence levels)
    to exclude when the genes are obtained from a \link{TranscriptDb} object.
  }
}

\value{
  A \link[Biostrings]{DNAStringSet} object containing one upstream sequence
  per gene (or per transcript if \code{genes} is a
  \link[GenomicRanges]{GenomicRanges} object containing transcript ranges).

  More precisely, if \code{genes} is a \link[GenomicRanges]{GenomicRanges}
  object, the returned object is \emph{parallel} to it, that is, the i-th
  element in the returned object is the upstream sequence corresponding to
  the i-th gene (or transcript) in \code{genes}. Also the names on the
  \link[GenomicRanges]{GenomicRanges} object are propagated to the returned
  object.

  If \code{genes} is a \link{TranscriptDb} object, the names on the returned
  object are the gene IDs found in the \link{TranscriptDb} object. To see the
  type of gene IDs (i.e. Entrez gene ID or Ensembl gene ID or ...), you can
  display \code{genes} with \code{show(genes)}.

  In addition, the returned object has the following metadata columns
  (accessible with \code{\link{mcols}}) that provide some information about
  the gene (or transcript) corresponding to each upstream sequence:
  \itemize{
    \item \code{gene_seqnames}: the chromosome name of the gene (or
          transcript);
    \item \code{gene_strand}: the strand of the gene (or transcript);
    \item \code{gene_TSS}: the transcription start site of the gene (or
          transcript).
  }
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item The \link[BSgenome]{BSgenome}, \link[rtracklayer]{TwoBitFile}, and
          \link[Rsamtools]{FaFile} classes, defined and documented
          in the \pkg{BSgenome}, \pkg{rtracklayer}, and \pkg{Rsamtools}
          packages, respectively.

    \item The \code{\link[BSgenome]{available.genomes}} function in the
          \pkg{BSgenome} package for checking avaibility of BSgenome
          data packages (and installing the desired one).
          
    \item The \link{TranscriptDb} class.

    \item The \code{\link{genes}} function for extracting gene ranges from
          a \link{TranscriptDb} object.

    \item The \link[GenomicRanges]{GenomicRanges} class defined and documented
          in the \pkg{GenomicRanges} package.

    \item The \link[Biostrings]{DNAStringSet} class defined and documented
          in the \pkg{Biostrings} package.

    \item The \code{\link[GenomeInfoDb]{seqinfo}} getter defined and documented
          in the \pkg{GenomeInfoDb} package.

    \item The \code{\link[Biostrings]{getSeq}} function for extracting
          subsequences from a sequence container.
  }
}

\examples{
## Load a genome:
library(BSgenome.Dmelanogaster.UCSC.dm3)
genome <- BSgenome.Dmelanogaster.UCSC.dm3
genome

## Use a TranscriptDb object:
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
txdb  # contains Ensembl gene IDs

## Because the chrU and chrUextra sequences are made of concatenated
## scaffolds (see http://genome.ucsc.edu/cgi-bin/hgGateway?db=dm3),
## extracting the upstream sequences for genes located on these
## scaffolds is not reliable. So we exclude them:
exclude <- c("chrU", "chrUextra")
up1000seqs <- extractUpstreamSeqs(genome, txdb, width=1000,
                                  exclude.seqlevels=exclude)
up1000seqs  # the names are Ensembl gene IDs
mcols(up1000seqs)

## Upstream sequences for genes close to the chromosome bounds can be
## shorter than 1000 (note that this does not happen for circular
## chromosomes like chrM):
table(width(up1000seqs))
mcols(up1000seqs)[width(up1000seqs) != 1000, ]
}

\keyword{manip}
