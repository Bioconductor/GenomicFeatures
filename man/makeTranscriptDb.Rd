\name{makeTranscriptDb}

\alias{GenomicFeatures-class}
\alias{class:GenomicFeatures}
\alias{GenomicFeatures}

\alias{TranscriptDb-class}
\alias{class:TranscriptDb}
\alias{TranscriptDb}

\alias{makeTranscriptDb}

\alias{show,TranscriptDb-method}
\alias{as.list,TranscriptDb-method}

\title{
  Making a TranscriptDb object from user supplied annotations
}
\description{
  \code{makeTranscriptDb} is a low-level constructor for making
  a TranscriptDb object from user supplied transcript annotations.
  See \code{?\link{makeTranscriptDbFromUCSC}} and
  \code{?\link{makeTranscriptDbFromBiomart}} for higher-level
  functions that feed data from the UCSC or BioMart sources
  to \code{makeTranscriptDb}.
}
\usage{
  makeTranscriptDb(transcripts, splicings, genes=NULL, metadata=NULL, ...)
}
\arguments{
  \item{transcripts}{data frame containing the genomic locations of
    a set of transcripts}
  \item{splicings}{data frame containing the exon and cds locations
    of a set of transcripts}
  \item{genes}{data frame containing the genes associated to a set
    of transcripts}
  \item{metadata}{2-column data frame containing meta information about
    this set of transcripts like species, organism, genome, UCSC table, etc...
    The names of the columns must be \code{"name"} and \code{"value"}
    and their type must be character.}
  \item{...}{ignored for now}
}
\details{
  The \code{transcripts}, \code{splicings} and \code{genes} arguments
  must be data frames that describe a set of transcripts and the genomic
  features related to them (exons, cds and genes at the moment).

  \code{transcripts} must have 1 row per transcript and the following
  columns:
  \itemize{
  \item \code{tx_id}: Transcript ID. Integer vector. No NAs. No duplicates. 
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
        Same type as \code{transcripts$tx_id} (integer vector). No NAs.
        All the values in this column must be present in
        \code{transcripts$tx_id}.
  \item \code{exon_rank}: The rank of the exon in the transcript.
        Integer vector with no NAs. (\code{tx_id}, \code{exon_rank})
        pairs must be unique.
  \item \code{exon_id}: [optional] Exon ID.
        Integer vector. No NAs.
  \item \code{exon_name}: [optional] Exon name.
        Character vector (or factor).
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
  \item \code{cds_id}: [optional] cds ID. Integer vector.
        If present then \code{cds_start} and \code{cds_end} must be too.
        NAs are allowed and must match NAs in \code{cds_start}
        and \code{cds_end}.
  \item \code{cds_name}: [optional] cds name. Character vector (or factor).
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

\value{A TranscriptDb object.}

\author{
  H. Pages
}

\seealso{
  \code{\link{makeTranscriptDbFromUCSC}},
  \code{\link{makeTranscriptDbFromBiomart}}
}

\examples{
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

txdb <- makeTranscriptDb(transcripts, splicings)

## Dump the entire TranscriptDb object into a list of data frames 'txdump'
## that can be used in 'do.call(makeTranscriptDb, txdump)' to make the db
## again with no loss of information.
## Note that the transcripts are dumped in a consistent order in all the
## elements of the result.
txdump <- as.list(txdb)
txdump
txdb1 <- do.call(makeTranscriptDb, txdump)
stopifnot(identical(as.list(txdb1), txdump))
}
