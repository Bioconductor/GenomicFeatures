\name{bindTranscripts}

\alias{bindTranscripts}
\alias{bindExons}
\alias{bindCDS}

\title{
  Binds genomic features from a TranscriptDb object to a GRanges object
}
\description{
  Functions to bind the features in a TranscriptDb object to a set of ranges.
}
\usage{
  bindTranscripts(txdb, ranges,
                  restrict=c("any", "start", "end", "within", "equal"),
                  columns=c("tx_id","tx_name"))
  bindExons(txdb, ranges, restrict=c("any", "start", "end", "within", "equal"),
            columns="exon_id")
  bindCDS(txdb, ranges, restrict=c("any", "start", "end", "within", "equal"),
         columns="cds_id")
}
\arguments{
  \item{txdb}{A \link{TranscriptDb} object.}
  \item{ranges}{A \link[GenomicRanges]{GRanges} object to restrict the output.}
  \item{restrict}{How to perform the interval overlap operations of the
    \code{ranges}. See the
    \code{\link[IRanges:IntervalTree-class]{findOverlaps}} manual page
    in the IRanges package for more information.}
  \item{columns}{columns to bind to \code{ranges}.
    For \code{bindTranscripts} can be NULL or a combination of "tx_id",
    "tx_name", "tx_chrom", "tx_ranges", and "tx_strand".
    For \code{bindExons} can be NULL or a combination of "exon_id",
    "exon_chrom", "exon_ranges", and "exon_strand".
    For \code{bindCDS} can be NULL or a combination of "cds_id",
    "cds_chrom", "cds_ranges", and "cds_strand".}
}
\details{
  These functions map specified \code{columns} that result from
  \code{\link{transcripts}}, \code{\link{exons}}, and \code{\link{cds}}
  calls to the input \code{ranges} using
  \code{\link[IRanges:IntervalTree-class]{findOverlaps}}.

  For \code{bindTranscripts}, bound \code{columns} consist of
  "tx_id" <CompressedIntegerList>,
  "tx_name" <CompressedCharacterList>,
  "tx_chrom" <CompressedCharacterList>,
  "tx_ranges" <CompressedIRangesList>, and
  "tx_strand" <CompressedCharacterList>.

  For \code{bindExons}, bound \code{columns} consist of
  "exon_id" <CompressedIntegerList>,
  "exon_chrom" <CompressedCharacterList>,
  "exon_ranges" <CompressedIRangesList>, and
  "exon_strand" <CompressedCharacterList>.

  For \code{bindCDS}, bound \code{columns} consist of
  "cds_id" <CompressedIntegerList>,
  "cds_chrom" <CompressedCharacterList>,
  "cds_ranges" <CompressedIRangesList>, and 
  "cds_strand" <CompressedCharacterList>.
}
\value{ a \code{\linkS4class{GRanges}} object }
\author{
  P. Aboyoun and M. Carlson
}
\seealso{
  \code{\linkS4class{TranscriptDb}}, \code{\link{transcripts}},
  \code{\link{transcriptsByOverlaps}}
}
\examples{
  txdb <- loadFeatures(system.file("extdata", "UCSC_knownGene_sample.sqlite",
                                   package="GenomicFeatures"))
  gr <- GRanges(seqnames = rep("chr1",2),
                       ranges = IRanges(start=c(500,10500),end=c(10000,30000)),
                       strand = rep("-",2))
  bindTranscripts(txdb, gr)
}
