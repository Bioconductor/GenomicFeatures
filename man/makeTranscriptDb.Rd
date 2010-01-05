\name{makeTranscriptDb}

\alias{makeTranscriptDb}
\alias{makeTranscriptDbFromUCSC}
\alias{makeTranscriptDbFromBiomart}

\alias{GenomicFeatures-class}
\alias{class:GenomicFeatures}
\alias{GenomicFeatures}

\alias{TranscriptDb-class}
\alias{class:TranscriptDb}
\alias{TranscriptDb}

\title{
  Creating TranscriptDb objects from publicly available sources.
}
\description{
  Flexible tools for creating a TranscriptDb object from
  a publicly available source (like the UCSC genome browser or BioMart),
  or from a local file.
}
\usage{
  makeTranscriptDb(geneId, txId, chrom, strand, txStart, txEnd,
                   cdsStart, cdsEnd, exonStart,exonEnd, exonRank, 
                   exonId)
  makeTranscriptDbFromUCSC(track=c("knownGene","refGene","ensGene"),
                           genome="hg18")
  makeTranscriptDbFromBiomart(biomart, dataset)
}
\arguments{
  \item{geneId}{vector of gene IDs}
  \item{txId}{vector of transcript IDs}
  \item{chrom}{vector of chromosomes}
  \item{strand}{vector of strand orientation}
  \item{txStart}{vector of transcription start sites}
  \item{txEnd}{vector of transcription end sites}
  \item{cdsStart}{vector of coding sequence start sites}
  \item{cdsEnd}{vector of coding sequence end sites}
  \item{exonStart}{vector of exon start sites}
  \item{exonEnd}{vector of exon end sites}
  \item{exonRank}{vector of exon ranks (position within a transcript)}
  \item{exonId}{vector of exon IDs}
  \item{track}{name of the UCSC track containing the transcript annotations
    to retrieve. One of "refGene", "knownGene", or "ensGene".}
  \item{genome}{genome abbreviation used by UCSC and obtained by
    \code{ucscGenomes()[,"db"]}. For example: hg18.}

  \item{biomart}{which mart from biomaRt}
  \item{dataset}{which dataset from biomaRt}
}
\details{
  All three of these functions perform similar tasks.

  The most general form is \code{makeTranscriptDb}, which simply takes
  vectors of equal length. The i-th element of the vectors represent data
  for a single exon. Transcript data must be repeated over all exons within
  a given trancsript.

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

txdb2 <- makeTranscriptDbFromBiomart("ensembl", "hsapiens_gene_ensembl")
}
}
