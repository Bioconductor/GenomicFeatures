\name{makeTranscriptDb}

\alias{GenomicFeatures-class}
\alias{class:GenomicFeatures}
\alias{GenomicFeatures}

\alias{TranscriptDb-class}
\alias{class:TranscriptDb}
\alias{TranscriptDb}

\alias{supportedUCSCtables}
\alias{makeTranscriptDbFromUCSC}
\alias{makeTranscriptDbFromBiomart}
\alias{makeTranscriptDb}

\alias{as.list,TranscriptDb-method}
\alias{as.data.frame,TranscriptDb-method}

\title{
  Creating TranscriptDb objects from publicly available sources.
}
\description{
  Flexible tools for creating a TranscriptDb object from
  a publicly available source (like the UCSC genome browser or BioMart),
  or from user supplied annotations.
}
\usage{
  supportedUCSCtables()
  makeTranscriptDbFromUCSC(genome="hg18", tablename="knownGene")
  makeTranscriptDbFromBiomart(biomart, dataset, ensembl_transcript_ids=NULL)
  makeTranscriptDb(transcripts, splicings, genes=NULL, ...)
}
\arguments{
  \item{genome}{genome abbreviation used by UCSC and obtained by
    \code{ucscGenomes()[ , "db"]}. For example: hg18.}
  \item{tablename}{name of the UCSC table containing the transcript
    annotations to retrieve. Use the \code{supportedUCSCtables} utility
    function to get the list of supported tables. Note that not all
    tables are available for all genomes.}
  \item{biomart}{which mart from BioMart}
  \item{dataset}{which dataset from BioMart}
  \item{ensembl_transcript_ids}{optionally, only retrieve transcript
    annotation data for the specified set of Ensembl transcript ids}
  \item{transcripts}{data frame containing the genomic locations of
    a set of transcripts}
  \item{splicings}{data frame containing the exon and cds locations
    of a set of transcripts}
  \item{genes}{data frame containing the genes associated to a set
    of transcripts}
  \item{...}{ignored for now}
}
\details{
  \code{makeTranscriptDbFromUCSC} and \code{makeTranscriptDbFromBiomart}
  are convenience functions that feed data from the UCSC or BioMart sources
  to the lower level \code{makeTranscriptDb} function.

  The \code{transcripts}, \code{splicings} and \code{genes} arguments
  of the \code{makeTranscriptDb} function must be data frames that
  describe a set of transcripts and the genomic features related to
  them (exons, cds and genes for now).

  \code{transcripts} must have 1 row per transcript and the following
  columns:
  \itemize{
  \item \code{tx_id}: Transcript ID. Character vector (or factor) or
        integer vector. No NAs. No duplicates. 
  \item \code{tx_name}: [optional] Transcript name. Character vector (or factor).
  \item \code{tx_chrom}: Transcript chromosome. Character vector (or factor)
        with no NAs.
  \item \code{tx_strand}: Transcript strand. Character vector (or factor)
        where each element is either \code{"+"} or \code{"-"}.
  \item \code{tx_start}, \code{tx_end}: Transcript start and end.
        Integer vectors with no NAs.
  }
  Other columns, if any, are ignored (with a warning).

  \code{splicings} must have N rows per transcript, where N is the nb
  of exons in the transcript. Each row describes an exon plus eventually
  the cds contained in this exon. Its columns must be:
  \itemize{
  \item \code{tx_id}: Foreign key that links each row in the \code{splicings}
        data frame to a unique row in the \code{transcripts} data frame.
        Note that more than 1 row in \code{splicings} can be linked to the
        same row in \code{transcripts} (many-to-one relationship).
        Same type as \code{transcripts$tx_id}. No NAs. All the values in
        this column must be present in \code{transcripts$tx_id}.
  \item \code{exon_rank}: The rank of the exon in the transcript.
        Integer vector with no NAs. (\code{tx_id}, \code{exon_rank})
        pairs must be unique.
  \item \code{exon_id}: [optional] Exon ID.
        Character vector (or factor) or integer vector. No NAs.
  \item \code{exon_chrom}: [optional] Exon chromosome.
        Character vector (or factor) with no NAs.
        If missing then \code{transcripts$tx_chrom} is used.
        If present then \code{exon_strand} must be present too.
  \item \code{exon_strand}: [optional] Exon strand.
        Character vector (or factor) with no NAs.
        If missing then \code{transcripts$tx_strand} is used
        and \code{exon_chrom} must be missing too.
  \item \code{exon_start}, \code{exon_end}: Exon start and end.
        Integer vectors with no NAs.
  \item \code{cds_id}: [optional] cds ID.
        Character vector (or factor) or integer vector.
        If present then \code{cds_start} and \code{cds_end} must be too.
        NAs are allowed and must match NAs in \code{cds_start}
        and \code{cds_end}.
  \item \code{cds_start}, \code{cds_end}: [optional] cds start and end.
        Integer vectors.
        If one of the 2 columns is missing then all \code{cds_*} columns
        must be missing.
        NAs are allowed and must occur at the same positions in
        \code{cds_start} and \code{cds_end}.
  }
  Other columns, if any, are ignored (with a warning).

  \code{genes} must have N rows per transcript, where N is the nb
  of genes linked to the transcript (N will be 1 most of the time).
  Its columns must be:
  \itemize{
  \item \code{tx_id}: [optional] \code{genes} must have either a
        \code{tx_id} or a \code{tx_name} column but not both.
        Like \code{splicings$tx_id}, this is a foreign key that
        links each row in the \code{genes} data frame to a unique
        row in the \code{transcripts} data frame.
  \item \code{tx_name}: [optional] 
        Can be used as an alternative to the \code{genes$tx_id}
        foreign key.
  \item \code{gene_id}: Gene ID. Character vector (or factor). No NAs.
  }
  Other columns, if any, are ignored (with a warning).
}

\value{All of these functions return a TranscriptDb object.}

\author{
  M. Carlson and H. Pages
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

## From user supplied annotations.
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

txdb0 <- makeTranscriptDb(transcripts, splicings)

## Dump the entire TranscriptDb object into a list of data frames 'txdump'
## that can be used in 'do.call(makeTranscriptDb, txdump)' to make the db
## again with no loss of information.
## Note that the transcripts are dumped in a consistent order in all the
## elements of the result.
txdump0 <- as.list(txdb0)
txdump0
txdb0b <- do.call(makeTranscriptDb, txdump0)
stopifnot(identical(as.list(txdb0b), txdump0))
}
