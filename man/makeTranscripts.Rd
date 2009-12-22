\name{makeTranscripts}
\alias{makeTranscripts}
\alias{UCSCTranscripts}
\alias{BMTranscripts}
\title{
  Methods to generate Transcripts objects from publicly available sources.
}
\description{
  Methods that will help to parse data and then assemble it into a
  coherent database so that it can be used as a Transcripts object.
}
\usage{
  makeTranscripts(geneIds, ids, chrom, strand, txStart, txEnd,
                  cdsStart, cdsEnd, exonStart,exonEnd, exonRank, 
                  exonId)
  UCSCTranscripts(type=c("refGene","knownGene","ensGene"), genome)
  BMTranscripts(biomart, dataset)
}
\arguments{
  \item{geneIds}{vector of GeneIds}
  \item{ids}{vector of transcript IDs}
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

  \item{type}{transcript table name for the UCSC track.
    One of "refGene", "knownGene", or "ensGene".}
  \item{genome}{genome abbreviation used by UCSC and obtained by
    \code{ucscGenomes()[,"db"]}. For example: hg18.}

  \item{biomart}{which mart from biomaRt}
  \item{dataset}{which dataset from biomaRt}
}
\details{
  All three of these functions perform similar tasks.  The most general
  form is makeTranscripts, which simply takes vectors of equal length
  and will then attempt to construct them into a transcript database.

  UCSCTranscripts and BMTranscripts are convenience functions for
  feeding data from these sources to the makeTranscripts function.
}
\value{ All of these functions return a Transcript Object. }
\author{
  M. Carlson
}

\examples{
\dontrun{
tx <- UCSCTranscripts()

tx <- BMTranscripts("ensembl", "hsapiens_gene_ensembl")
}
}
