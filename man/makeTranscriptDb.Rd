\name{makeTranscriptDb}

\alias{GenomicFeatures-class}
\alias{class:GenomicFeatures}
\alias{GenomicFeatures}

\alias{TranscriptDb-class}
\alias{class:TranscriptDb}
\alias{TranscriptDb}

\alias{makeTranscriptDb}
\alias{makeTranscriptDbFromUCSC}
\alias{makeTranscriptDbFromBiomart}

\alias{as.data.frame,TranscriptDb-method}

\title{
  Creating TranscriptDb objects from publicly available sources.
}
\description{
  Flexible tools for creating a TranscriptDb object from
  a publicly available source (like the UCSC genome browser or BioMart),
  or from a local file.
}
\usage{
  makeTranscriptDb(geneId,
                   txId, txName=NULL, txChrom, txStrand, txStart, txEnd,
                   cdsStart, cdsEnd, exonStart,exonEnd, exonRank, 
                   exonId=NULL)
  makeTranscriptDbFromUCSC(genome="hg18",
                           track=c("knownGene","refGene","ensGene"))
  makeTranscriptDbFromBiomart(biomart, dataset, ensembl_transcript_ids=NULL)
}
\arguments{
  \item{geneId}{vector of gene IDs}
  \item{txId}{vector of transcript IDs}
  \item{txName}{vector of transcript names (not necessarily unique)}
  \item{txChrom}{vector of chromosomes}
  \item{txStrand}{vector of strand orientation}
  \item{txStart}{vector of transcription start sites}
  \item{txEnd}{vector of transcription end sites}
  \item{cdsStart}{vector of coding sequence start sites}
  \item{cdsEnd}{vector of coding sequence end sites}
  \item{exonStart}{vector of exon start sites}
  \item{exonEnd}{vector of exon end sites}
  \item{exonRank}{vector of exon ranks (position within a transcript)}
  \item{exonId}{vector of exon IDs}
  \item{genome}{genome abbreviation used by UCSC and obtained by
    \code{ucscGenomes()[,"db"]}. For example: hg18.}
  \item{track}{name of the UCSC track containing the transcript annotations
    to retrieve. One of "refGene", "knownGene", or "ensGene".}

  \item{biomart}{which mart from BioMart}
  \item{dataset}{which dataset from BioMart}
  \item{ensembl_transcript_ids}{optionally, only retrieve transcript
    annotation data for the specified set of Ensembl transcript ids}
}
\details{
  All three of these functions perform similar tasks.

  The most general form is \code{makeTranscriptDb}, which simply takes
  vectors of equal length. The i-th element of the vectors represent data
  for a single exon. Transcript data must be repeated over all exons within
  a given transcript.

  \code{makeTranscriptDbFromUCSC} and \code{makeTranscriptDbFromBiomart}
  are convenience functions for feeding data from these sources (UCSC or
  BioMart) to the \code{makeTranscriptDb} function.
}
\value{ All of these functions return a TranscriptDb object. }
\author{
  M. Carlson
}

\examples{
\dontrun{
txdb1 <- makeTranscriptDbFromUCSC()

txdb2 <- makeTranscriptDbFromBiomart()
}
## Make a tiny TranscriptDb object:
some_ids <- c(
    "ENST00000400839",
    "ENST00000400840",
    "ENST00000478783",
    "ENST00000435657",
    "ENST00000268655",
    "ENST00000313243",
    "ENST00000341724"
)
tiny_txdb <- makeTranscriptDbFromBiomart(ensembl_transcript_ids=some_ids)

### Low-level utility.
transcripts <- data.frame(
                   tx_id=1:3,
                   tx_chrom="chr1",
                   tx_strand=c("-", "+", "+"),
                   tx_start=c(1, 2001, 2001),
                   tx_end=c(999, 2199, 2199))
splicings <-  data.frame(
                   tx_id=c(1L, 2L, 2L, 2L, 3L, 3L),
                   exon_rank=c(1, 1, 2, 3, 1, 2),
                   exon_start=c(1, 2001, 2101, 2131, 2001, 2131),
                   exon_end=c(999, 2085, 2144, 2199, 2085, 2199),
                   cds_start=c(1, 2022, 2101, 2131, NA, NA),
                   cds_end=c(999, 2085, 2144, 2193, NA, NA))

txdb <- GenomicFeatures:::.makeTranscriptDb(transcripts, splicings)
}
