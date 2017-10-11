\name{exonicParts}

\alias{exonicParts}
\alias{intronicParts}

\title{
  Extract non-overlapping exonic or intronic parts from a TxDb-like object
}

\description{
  \code{exonicParts} and \code{intronicParts} extract the non-overlapping
  (a.k.a. disjoint) exonic or intronic parts from a \link{TxDb}-like object.
}

\usage{
exonicParts(txdb, linked.to.single.gene.only=FALSE)
intronicParts(txdb, linked.to.single.gene.only=FALSE)
}

\arguments{
  \item{txdb}{
    A \link{TxDb} object, or any \link{TxDb}-like object that supports the
    \code{\link{transcripts}()} and \code{\link{exonsBy}()} extractors
    (e.g. an \link[ensembldb]{EnsDb} object).
  }
  \item{linked.to.single.gene.only}{
    \code{TRUE} or \code{FALSE}.

    If \code{FALSE} (the default), then the disjoint parts are obtained
    by calling \code{\link{disjoin()}} on all the exons (or introns) in
    \code{txdb}, including on exons (or introns) not linked to a gene or
    linked to more than one gene.

    If \code{TRUE}, then the disjoint parts are obtained in 2 steps:
    \enumerate{
      \item call \code{\link{disjoin()}} on the exons (or introns)
            linked to \emph{at least one gene},

      \item then drop the parts linked to more than one gene from
            the set of exonic (or intronic) parts obtained previously.
    }
  }
}

\value{
  \code{exonicParts} returns a disjoint and strictly sorted
  \link[GenomicRanges]{GRanges} object with 1 range per exonic part
  and with metadata columns \code{tx_id}, \code{tx_name}, \code{gene_id},
  \code{exon_id}, \code{exon_name}, and \code{exon_rank}.

  \code{intronicParts} returns a disjoint and strictly sorted
  \link[GenomicRanges]{GRanges} object with 1 range per intronic part
  and with metadata columns \code{tx_id}, \code{tx_name}, and \code{gene_id}.
}

\note{
  \code{exonicParts} is a replacement for \code{\link{disjointExons}} with
  the following differences/improvements:
  \itemize{
    \item Argument \code{linked.to.single.gene.only} in \code{exonicParts}
          replaces argument \code{aggregateGenes} in \code{disjointExons},
          but has opposite meaning i.e.
          \code{exonicParts(txdb, linked.to.single.gene.only=TRUE)}
          returns the same exonic parts as
          \code{disjointExons(txdb, aggregateGenes=FALSE)}.

    \item Unlike \code{disjointExons(txdb, aggregateGenes=TRUE)},
          \code{exonicParts(txdb, linked.to.single.gene.only=FALSE)} does
          NOT discard exon parts that are not linked to a gene.

    \item \code{exonicParts} is almost 2x more efficient than
          \code{disjointExons}.

    \item \code{exonicParts} works out-of-the-box on any \link{TxDb}-like
          object that supports the \code{\link{transcripts}()} and
          \code{\link{exonsBy}()} extractors (e.g. on an
          \link[ensembldb]{EnsDb} object).
  }
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \code{\link{transcripts}}, \code{\link{transcriptsBy}},
          and \code{\link{transcriptsByOverlaps}}, for extracting
          genomic feature locations from a \link{TxDb}-like object.

    \item \code{\link{transcriptLengths}} for extracting the transcript
          lengths (and other metrics) from a \link{TxDb} object.

    \item \code{\link{extractTranscriptSeqs}} for extracting transcript
          (or CDS) sequences from chromosome sequences.

    \item \code{\link{coverageByTranscript}} for computing coverage by
          transcript (or CDS) of a set of ranges.

    \item The \link{TxDb} class.
  }
}

\examples{
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## ---------------------------------------------------------------------
## exonicParts()
## ---------------------------------------------------------------------

exonic_parts1 <- exonicParts(txdb)
exonic_parts1

## Mapping from exonic parts to genes is many-to-many:
mcols(exonic_parts1)$gene_id
table(lengths(mcols(exonic_parts1)$gene_id))
## A Human exonic part can be linked to 0 to 22 known genes!

exonic_parts2 <- exonicParts(txdb, linked.to.single.gene.only=TRUE)
exonic_parts2

## Mapping from exonic parts to genes now is many-to-one:
class(mcols(exonic_parts2)$gene_id)

## Sanity checks:
stopifnot(isDisjoint(exonic_parts1), isStrictlySorted(exonic_parts1))
stopifnot(isDisjoint(exonic_parts2), isStrictlySorted(exonic_parts2))
stopifnot(all(exonic_parts2 \%within\% reduce(exonic_parts1)))
stopifnot(identical(
    lengths(mcols(exonic_parts1)$gene_id) == 1L,
    exonic_parts1 \%within\% exonic_parts2
))

## ---------------------------------------------------------------------
## intronicParts()
## ---------------------------------------------------------------------

intronic_parts1 <- intronicParts(txdb)
intronic_parts1

## Mapping from intronic parts to genes is many-to-many:
mcols(intronic_parts1)$gene_id
table(lengths(mcols(intronic_parts1)$gene_id))
## A Human intronic part can be linked to 0 to 22 known genes!

intronic_parts2 <- intronicParts(txdb, linked.to.single.gene.only=TRUE)
intronic_parts2

## Mapping from intronic parts to genes now is many-to-one:
class(mcols(intronic_parts2)$gene_id)

## Sanity checks:
stopifnot(isDisjoint(intronic_parts1), isStrictlySorted(intronic_parts1))
stopifnot(isDisjoint(intronic_parts2), isStrictlySorted(intronic_parts2))
stopifnot(all(intronic_parts2 \%within\% reduce(intronic_parts1)))
stopifnot(identical(
    lengths(mcols(intronic_parts1)$gene_id) == 1L,
    intronic_parts1 \%within\% intronic_parts2
))
}

\keyword{manip}
