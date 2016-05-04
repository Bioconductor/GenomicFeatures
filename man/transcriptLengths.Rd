\name{transcriptLengths}

\alias{transcriptLengths}


\title{Extract the transcript lengths from a TxDb object}

\description{
  The \code{transcriptLengths} function extracts the transcript lengths from
  a \link{TxDb} object. It also returns the CDS and UTR lengths for each
  transcript if the user requests them.
}

\usage{
transcriptLengths(txdb, with.cds_len=FALSE,
                  with.utr5_len=FALSE, with.utr3_len=FALSE, ...)
}

\arguments{
  \item{txdb}{
    A \link{TxDb} object.
  }
  \item{with.cds_len, with.utr5_len, with.utr3_len}{
    \code{TRUE} or \code{FALSE}. Whether or not to also extract and return
     the CDS, 5' UTR, and 3' UTR lengths for each transcript.
  }
  \item{\dots}{
    Additional arguments used by \code{transcripts} and other accessor 
    functions.
  }
}

\details{
  All the lengths are counted in number of nucleotides.

  The length of a processed transcript is just the sum of the lengths of its
  exons. This should not be confounded with the length of the stretch of DNA
  transcribed into RNA (a.k.a. transcription unit), which can be obtained
  with \code{width(transcripts(txdb))}.
}

\value{
  A data frame with 1 row per transcript. The rows are guaranteed to be in
  the same order as the elements of the \link[GenomicRanges]{GRanges} object
  returned by \code{\link{transcripts}(txdb)}.
  The data frame has between 5 and 8 columns, depending on what the user
  requested via the \code{with.cds_len}, \code{with.utr5_len}, and
  \code{with.utr3_len} arguments.

  The first 3 columns are the same as the metadata columns of the object
  returned by
\preformatted{  transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
}
  that is:
  \itemize{
    \item \code{tx_id}: The internal transcript ID. This ID is unique within
          the scope of the \link{TxDb} object. It is not an official or public
          ID (like an Ensembl or FlyBase ID) or an Accession number, so it
          cannot be used to lookup the transcript in public data bases or in
          other \link{TxDb} objects. Furthermore, this ID could change when
          re-running the code that was used to make the \link{TxDb} object.
    \item \code{tx_name}: An official/public transcript name or ID that can
          be used to lookup the transcript in public data bases or in other
          \link{TxDb} objects. This column is not guaranteed to contain unique
          values and it can contain NAs.
    \item \code{gene_id}: The official/public ID of the gene that the
          transcript belongs to. Can be NA if the gene is unknown or if the
          transcript is not considered to belong to a gene.
  }

  The other columns are quantitative:
  \itemize{
    \item \code{nexon}: The number of exons in the transcript.
    \item \code{tx_len}: The length of the processed transcript.
    \item \code{cds_len}: [optional] The length of the CDS region of the
          processed transcript.
    \item \code{utr5_len}: [optional] The length of the 5' UTR region of
          the processed transcript.
    \item \code{utr3_len}: [optional] The length of the 3' UTR region of
          the processed transcript.
  }
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \code{\link{transcripts}}, \code{\link{transcriptsBy}},
          and \code{\link{transcriptsByOverlaps}},
          for how to extract the genomic locations of features from a
          \link{TxDb} object.

    \item \code{\link{extractTranscriptSeqs}} for extracting transcript
          (or CDS) sequences from chromosome sequences.

    \item \code{\link{coverageByTranscript}} for computing coverage by
          transcript (or CDS) of a set of ranges.

    \item \code{\link{makeTxDbFromUCSC}} and
          \code{\link{makeTxDbFromBiomart}} for convenient ways to
          make \link{TxDb} objects from UCSC or BioMart online resources.

    \item The \link{TxDb} class.
  }
}

\examples{
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
dm3_txlens <- transcriptLengths(txdb)
head(dm3_txlens)

dm3_txlens <- transcriptLengths(txdb, with.cds_len=TRUE,
                                      with.utr5_len=TRUE,
                                      with.utr3_len=TRUE)
head(dm3_txlens)

## When cds_len is 0 (non-coding transcript), utr5_len and utr3_len
## must also be 0:
non_coding <- dm3_txlens[dm3_txlens$cds_len == 0, ]
stopifnot(all(non_coding[6:8] == 0))

## When cds_len is not 0 (coding transcript), cds_len + utr5_len +
## utr3_len must be equal to tx_len:
coding <- dm3_txlens[dm3_txlens$cds_len != 0, ]
stopifnot(all(rowSums(coding[6:8]) == coding[[5]]))

## A sanity check:
stopifnot(identical(dm3_txlens$tx_id, mcols(transcripts(txdb))$tx_id))
}

\keyword{manip}
