\name{disjointExons}

\alias{disjointExons}
\alias{disjointExons,TxDb-method}

\title{
  Extract non-overlapping exon parts from an object
}

\description{
  \code{disjointExons} extracts the non-overlapping exon parts from a
  \link{TxDb} object or any other supported object.
}

\usage{
disjointExons(x, ...)

\S4method{disjointExons}{TxDb}(x, aggregateGenes=FALSE, 
                 includeTranscripts=TRUE, ...) 
}

\arguments{ 
  \item{x}{
    A \link{TxDb} object or any other supported object.
  }
  \item{...}{
    Arguments to be passed to methods.
  }
  \item{aggregateGenes}{
    For \code{disjointExons} : A \code{logical}. When \code{FALSE}
    (default) exon fragments that overlap multiple genes are dropped.
    When \code{TRUE}, all fragments are kept and the \code{gene_id}
    metadata column includes all gene ids that overlap the exon fragment.
  }
  \item{includeTranscripts}{
    For \code{disjointExons} : A \code{logical}. When \code{TRUE}
    (default) a \code{tx_name} metadata column is included that 
    lists all transcript names that overlap the exon fragment. 
  }
}

\details{
  \code{disjointExons} creates a \link[GenomicRanges]{GRanges} of
  non-overlapping exon parts with metadata columns of gene_id and exonic_part.
  Exon parts that overlap more than 1 gene can be dropped with 
  \code{aggregateGenes=FALSE}. When \code{includeTranscripts=TRUE} 
  a \code{tx_name} metadata column is included that lists all 
  transcript names that overlap the exon fragment. This function 
  replaces \code{prepareAnnotationForDEXSeq} in the \pkg{DEXSeq} 
  package.
}

\value{
  A \link[GenomicRanges]{GRanges} object.
}

\author{
  \code{disjointExons} was originally implemented by Mike Love and Alejandro
  Reyes and then moved (and adapted) to \pkg{GenomicFeatures} by Valerie
  Obenchain.
}

\seealso{
  \itemize{
    \item \code{\link{transcripts}}, \code{\link{transcriptsBy}}, and
          \code{\link{transcriptsByOverlaps}} for the core genomic features
          extractors.
    \item The \link{TxDb} class.
  }
}

\examples{
## TODO
}

\keyword{methods}
