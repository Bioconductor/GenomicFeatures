\name{transcriptLocs2refLocs}

\alias{transcriptWidths}
\alias{transcriptLocs2refLocs}

\title{Converting transcript-based locations into reference-based locations}

\description{
  \code{transcriptLocs2refLocs} converts transcript-based
  locations into reference-based (aka chromosome-based or genomic)
  locations.

  \code{transcriptWidths} computes the lengths of the transcripts
  (called the "widths" in this context) based on the boundaries
  of their exons.
}

\usage{
transcriptLocs2refLocs(tlocs,
        exonStarts=list(), exonEnds=list(), strand=character(0),
        decreasing.rank.on.minus.strand=FALSE, error.if.out.of.bounds=TRUE)

transcriptWidths(exonStarts=list(), exonEnds=list())
}

\arguments{
  \item{tlocs}{
    A list of integer vectors of the same length as \code{exonStarts}
    and \code{exonEnds}. Each element in \code{tlocs} must contain
    transcript-based locations.
  }
  \item{exonStarts, exonEnds}{
    The starts and ends of the exons, respectively.

    Each argument can be a list of integer vectors,
    an \link[IRanges]{IntegerList} object,
    or a character vector where each element is a
    comma-separated list of integers.
    In addition, the lists represented by \code{exonStarts}
    and \code{exonEnds} must have the same shape i.e.
    have the same lengths and have elements of the same lengths.
    The length of \code{exonStarts} and \code{exonEnds}
    is the number of transcripts.
  }
  \item{strand}{
    A character vector of the same length as \code{exonStarts} and
    \code{exonEnds} specifying the strand (\code{"+"} or \code{"-"})
    from which the transcript is coming.
  }
  \item{decreasing.rank.on.minus.strand}{
    \code{TRUE} or \code{FALSE}.
    Describes the order of exons in transcripts located on the minus strand:
    are they ordered by increasing (default) or decreasing rank?
  }
  \item{error.if.out.of.bounds}{
    \code{TRUE} or \code{FALSE}.
    Controls how out of bound \code{tlocs} are handled: an error is thrown
    (default) or \code{NA} is returned.
  }
}

\value{
  For \code{transcriptLocs2refLocs}: A list of integer vectors of the same
  shape as \code{tlocs}.

  For \code{transcriptWidths}: An integer vector with one element per
  transcript.
}

\author{Hervé Pagès}

\seealso{
  \itemize{
    \item \code{\link{extractTranscriptSeqs}} for extracting transcript
          (or CDS) sequences from chromosomes.

    \item \code{\link{coverageByTranscript}} for computing coverage by
          transcript (or CDS) of a set of ranges.
  }
}

\examples{
## ---------------------------------------------------------------------
## WITH A SMALL SET OF HUMAN TRANSCRIPTS
## ---------------------------------------------------------------------
txdb_file <- system.file("extdata", "hg19_knownGene_sample.sqlite",
                         package="GenomicFeatures")
txdb <- loadDb(txdb_file)
ex_by_tx <- exonsBy(txdb, by="tx", use.names=TRUE)
genome <- BSgenome::getBSgenome("hg19")  # load the hg19 genome
tx_seqs <- extractTranscriptSeqs(genome, ex_by_tx)

## Get the reference-based locations of the first 4 (5' end)
## and last 4 (3' end) nucleotides in each transcript:
tlocs <- lapply(width(tx_seqs), function(w) c(1:4, (w-3):w))
tx_strand <- sapply(strand(ex_by_tx), runValue)

## Note that, because of how we made them, 'tlocs', 'start(ex_by_tx)',
## 'end(ex_by_tx)' and 'tx_strand' are "parallel" objects i.e. they
## have the same length, and, for any valid positional index, elements
## at this position are corresponding to each other. This is how
## transcriptLocs2refLocs() expects them to be!
rlocs <- transcriptLocs2refLocs(tlocs,
             start(ex_by_tx), end(ex_by_tx),
             tx_strand, decreasing.rank.on.minus.strand=TRUE)

## ---------------------------------------------------------------------
## WITH TWO WORM TRANSCRIPTS: ZC101.3.1 AND F37B1.1.1
## ---------------------------------------------------------------------
library(TxDb.Celegans.UCSC.ce11.ensGene)
txdb <- TxDb.Celegans.UCSC.ce11.ensGene
my_tx_names <- c("ZC101.3.1", "F37B1.1.1")
## Both transcripts are on chromosome II, the first one on its positive
## strand and the second one on its negative strand:
my_tx <- transcripts(txdb, filter=list(tx_name=my_tx_names))
my_tx

## Using transcripts stored in a GRangesList object:
ex_by_tx <- exonsBy(txdb, use.names=TRUE)[my_tx_names]
genome <- getBSgenome("ce11")  # load the ce11 genome
tx_seqs <- extractTranscriptSeqs(genome, ex_by_tx)
tx_seqs

## Since the 2 transcripts are on the same chromosome, an alternative
## is to store them in an IRangesList object and use that object with
## extractTranscriptSeqs():
ex_by_tx2 <- ranges(ex_by_tx)
tx_seqs2 <- extractTranscriptSeqs(genome$chrII, ex_by_tx2,
                                  strand=strand(my_tx))
stopifnot(identical(as.character(tx_seqs), as.character(tx_seqs2)))

## Store exon starts and ends in two IntegerList objects for use with
## transcriptWidths() and transcriptLocs2refLocs():
exon_starts <- start(ex_by_tx)
exon_ends <- end(ex_by_tx)

## Same as 'width(tx_seqs)':
transcriptWidths(exonStarts=exon_starts, exonEnds=exon_ends)

transcriptLocs2refLocs(list(c(1:2, 202:205, 1687:1688),
                            c(1:2, 193:196, 721:722)),
                       exonStarts=exon_starts,
                       exonEnds=exon_ends,
                       strand=c("+","-"))

## A sanity check:
ref_locs <- transcriptLocs2refLocs(list(1:1688, 1:722),
                                   exonStarts=exon_starts,
                                   exonEnds=exon_ends,
                                   strand=c("+","-"))
stopifnot(genome$chrII[ref_locs[[1]]] == tx_seqs[[1]])
stopifnot(complement(genome$chrII)[ref_locs[[2]]] == tx_seqs[[2]])
}

\keyword{manip}
