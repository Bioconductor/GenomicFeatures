\name{getPromoterSeq}

\alias{getPromoterSeq}
\alias{getTerminatorSeq}
\alias{getPromoterSeq,GRanges-method}
\alias{getTerminatorSeq,GRanges-method}
\alias{getPromoterSeq,GRangesList-method}
\alias{getTerminatorSeq,GRangesList-method}

\title{Get gene promoter or terminator sequences}

\description{
  Extract promoter or terminator sequences for the genes or transcripts
  specified in the query (a\link{GRanges} or \link{GRangesList} object)
  from a \link[BSgenome]{BSgenome} or \link[Rsamtools]{FaFile} object.
}

\usage{
\S4method{getPromoterSeq}{GRanges}(query, subject, upstream=2000, downstream=200)
\S4method{getTerminatorSeq}{GRanges}(query, subject, upstream=2000, downstream=200)

\S4method{getPromoterSeq}{GRangesList}(query, subject, upstream=2000, downstream=200)
\S4method{getTerminatorSeq}{GRangesList}(query, subject, upstream=2000, downstream=200)
}

\arguments{
  \item{query}{A \link[GenomicRanges]{GRanges} or
    \link[GenomicRanges]{GRangesList} object containing genes grouped by
    transcript.
  }
  \item{subject}{A \link[BSgenome]{BSgenome} or \link[Rsamtools]{FaFile} object from which
    the sequences will be taken.}
  \item{upstream}{The number of DNA bases to include upstream of the TSS (transcription start site)}
  \item{downstream}{The number of DNA bases to include downstream of the TSS (transcription start site)}
}

\details{
  \code{getPromoterSeq} and \code{getTerminatorSeq} are generic functions
  dispatching on query, which is either a GRanges or a GRangesList.
  They are convenience wrappers for the \code{promoters}, \code{terminators},
  and \code{getSeq} functions.
  The purpose is to allow sequence extraction from either a
  \link[BSgenome]{BSgenome} or \link[Rsamtools]{FaFile} object.

  Default values for \code{upstream} and \code{downstream} were chosen based
  on our current understanding of gene regulation. On average, promoter
  regions in the mammalian genome are 5000 bp upstream and downstream of the
  transcription start site.
}

\value{
  A \link[Biostrings]{DNAStringSet} or
  \link[Biostrings]{DNAStringSetList} instance corresponding to the
  GRanges or GRangesList supplied in the query.
}

\author{Paul Shannon}

\seealso{
  \itemize{
    \item The \code{\link[GenomicRanges]{promoters}} man page in the
          \pkg{GenomicRanges} package for the \code{promoters()} and
          \code{terminators()} methods for \link[GenomicRanges]{GenomicRanges}
          objects.

    \item \code{\link[Biostrings]{getSeq}} in the \pkg{Biostrings}
          package for extracting a set of sequences from a sequence
          container like a \link[BSgenome]{BSgenome} or
          \link[Rsamtools]{FaFile} object.
  }
}

\examples{
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)


## A GRangesList object describing all the known Human transcripts grouped
## by gene:
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tx_by_gene <- transcriptsBy(txdb, by="gene")

e2f3 <- "1871"  # entrez geneID for a cell cycle control transcription
                # factor, chr6 on the plus strand

## A GRanges object describing the three transcripts for gene 1871:
e2f3_tx <- tx_by_gene[[e2f3]]

## Promoter sequences for gene 1871:
e2f3_promoter_seqs <- getPromoterSeq(e2f3_tx, Hsapiens,
                                     upstream=40, downstream=15)
e2f3_promoter_seqs

mcols(e2f3_promoter_seqs)

## Terminator sequences for gene 1871:
e2f3_terminator_seqs <- getTerminatorSeq(e2f3_tx, Hsapiens,
                                         upstream=25, downstream=10)

e2f3_terminator_seqs

mcols(e2f3_terminator_seqs)  # same as 'mcols(e2f3_promoter_seqs)'

## All Human promoter sequences grouped by gene:
getPromoterSeq(tx_by_gene, Hsapiens, upstream=6, downstream=4)
}

\keyword{methods}
\keyword{manip}
