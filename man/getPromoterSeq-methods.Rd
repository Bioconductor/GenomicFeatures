\name{getPromoterSeq}
\alias{getPromoterSeq}
\alias{getPromoterSeq,GRangesList-method}
\alias{getPromoterSeq,GRanges-method}

\title{Get gene promoter sequences}

\description{
  Extract sequences for the genes or transcripts specified in the query
  (a\link{GRanges} or \link{GRangesList} object)
  from a \link[BSgenome]{BSgenome} object or an \link[Rsamtools]{FaFile}.
}

\usage{
  \S4method{getPromoterSeq}{GRangesList}(query, subject, upstream=2000, downstream=200, ...)
  \S4method{getPromoterSeq}{GRangesList}(query, subject, upstream=2000, downstream=200, ...)
  \S4method{getPromoterSeq}{GRanges}(query, subject, upstream=2000, downstream=200, ...)
}

\arguments{
  \item{query}{A \link[GenomicRanges]{GRanges} or \link[GenomicRanges]{GRangesList} object containing
  genes grouped by transcript.
  }
  \item{subject}{A \link[BSgenome]{BSgenome} object or a \link[Rsamtools]{FaFile} from which
    the sequences will be taken.}
  \item{upstream}{The number of DNA bases to include upstream of the TSS (transcription start site)}
  \item{downstream}{The number of DNA bases to include downstream of the TSS (transcription start site)}
  \item{\dots}{Additional arguments}
}

\details{
  \code{getPromoterSeq} is an overloaded method dispatching on query,
  which is either a GRanges or a GRangesList.   It is a wrapper for the
  \code{promoters} and \code{getSeq} functions. The 
  purpose is to allow sequence extraction from either a
  \link[BSgenome]{BSgenome} or \link[Rsamtools]{FaFile}.

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
  \link[IRanges]{intra-range-methods}        ## promoters methods for Ranges objects
  \link[GenomicRanges]{intra-range-methods}  ## promoters methods for GRanges objects
  \link[Biostrings]{getSeq}
}

\examples{
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg19)
 
  e2f3 <- "1871"  # entrez geneID for a cell cycle control transcription
                  # factor, chr6 on the plus strand

  transcriptCoordsByGene.GRangesList <-
     transcriptsBy (TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene") [e2f3]
    # a GrangesList of length one, describing three transcripts

  promoter.seqs <- getPromoterSeq (transcriptCoordsByGene.GRangesList,
                                   Hsapiens, upstream=10, downstream=0)
    # DNAStringSetList of length 1
    # [["1871"]] GCTTCCTGGA GCTTCCTGGA CGGAGCCAGG

}

\keyword{methods}
\keyword{manip}
