\name{extractTranscriptsFromGenome}

\alias{transcriptWidths}
\alias{transcriptLocs2refLocs}
\alias{extractTranscripts}
\alias{extractTranscriptsFromGenome}
\alias{sortExonsByRank}

\title{Tools for extracting transcript sequences}

\description{
  \code{extractTranscriptsFromGenome} extracts the transcript
  sequences from a BSgenome data package using the transcript
  information (exon boundaries) stored in a \link{TranscriptDb}
  or \link[GenomicRanges]{GRangesList} object.

  \code{extractTranscripts} extracts a set of transcripts
  from a single DNA sequence.

  Related utilities:

  \code{transcriptWidths} to get the lengths of the transcripts
  (called the "widths" in this context) based on the boundaries
  of their exons.

  \code{transcriptLocs2refLocs} converts transcript-based
  locations into reference-based (aka chromosome-based or genomic)
  locations.

  \code{sortExonsByRank} orders (or reorders) by rank the exons
  stored in a \link[GenomicRanges]{GRangesList} object containing
  exons grouped by transcript.
}

\usage{
extractTranscriptsFromGenome(genome, txdb,
                             decreasing.rank.on.minus.strand=FALSE,
                             use.names=TRUE)

extractTranscripts(x,
        exonStarts=list(), exonEnds=list(), strand=character(0),
        decreasing.rank.on.minus.strand=FALSE)

## Related utilities:

transcriptWidths(exonStarts=list(), exonEnds=list())

transcriptLocs2refLocs(tlocs,
        exonStarts=list(), exonEnds=list(), strand=character(0),
        decreasing.rank.on.minus.strand=FALSE)

sortExonsByRank(x, decreasing.rank.on.minus.strand=FALSE)
}

\arguments{
  \item{genome}{
    A \link[BSgenome]{BSgenome} object.
    See the \code{\link[BSgenome]{available.genomes}} function
    in the BSgenome package for how to install a genome.
  }
  \item{txdb}{
    A \link{TranscriptDb} object or a \link[GenomicRanges]{GRangesList}
    object.
  }
  \item{decreasing.rank.on.minus.strand}{
    \code{TRUE} or \code{FALSE}.
    Describes the order of exons in transcripts located on the minus strand:
    are they ordered by increasing (default) or decreasing rank?
    For all the functions described in this man page (except
    \code{sortExonsByRank}), this argument describes the input.
    For \code{sortExonsByRank}, it describes how exons should be ordered in
    the output.
  }
  \item{use.names}{
    \code{TRUE} or \code{FALSE}. Ignored if \code{txdb} is not
    a \link{TranscriptDb} object.
    If \code{TRUE} (the default), the returned sequences are named
    with the transcript names.
    If \code{FALSE}, they are named with the transcript internal ids.
    Note that, unlike the transcript internal ids, the transcript names
    are not guaranteed to be unique or even defined (they could be all
    \code{NA}s). A warning is issued when this happens.
  }
  \item{x}{
    A \link[Biostrings]{DNAString} or \link[Biostrings]{MaskedDNAString}
    object for \code{extractTranscripts}.

    A \link[GenomicRanges]{GRangesList} object for \code{sortExonsByRank},
    typically coming from \code{exonsBy(... , by="tx")}.
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
  \item{tlocs}{
    A list of integer vectors of the same length as \code{exonStarts}
    and \code{exonEnds}. Each element in \code{tlocs} must contain
    transcript-based locations.
  }
}

\value{
  For \code{extractTranscriptsFromGenome}: A named
  \link[Biostrings]{DNAStringSet} object with one element per transcript.
  When \code{txdb} is a \link[GenomicRanges]{GRangesList} object, elements
  in the output align with elements in the input (\code{txdb}), and they
  have the same names.

  For \code{extractTranscripts}: A \link[Biostrings]{DNAStringSet} object
  with one element per transcript.

  For \code{transcriptWidths}: An integer vector with one element per
  transcript.

  For \code{transcriptLocs2refLocs}: A list of integer vectors of the same
  shape as \code{tlocs}.

  For \code{sortExonsByRank}: A \link[GenomicRanges]{GRangesList} object
  with one top-level element per transcript. More precisely, the returned
  object has the same "shape" (i.e. same length and same number of elements
  per top-level element) as the input \link[GenomicRanges]{GRangesList}
  object \code{x}.
}

\author{
  H. Pages
}

\seealso{
  \code{\link[BSgenome]{available.genomes}},
  \link{TranscriptDb-class},
  \code{\link{exonsBy}},
  \link[GenomicRanges]{GRangesList-class},
  \link[Biostrings]{DNAStringSet-class},
  \code{\link[Biostrings]{translate}}
}

\examples{
library(BSgenome.Hsapiens.UCSC.hg18)  # load the genome

## ---------------------------------------------------------------------
## A. USING extractTranscriptsFromGenome() WITH A TranscriptDb OBJECT
## ---------------------------------------------------------------------
txdb_file <- system.file("extdata", "UCSC_knownGene_sample.sqlite",
                         package="GenomicFeatures")
txdb <- loadFeatures(txdb_file)
tx_seqs1 <- extractTranscriptsFromGenome(Hsapiens, txdb)
tx_seqs1

## ---------------------------------------------------------------------
## B. USING extractTranscriptsFromGenome() WITH A GRangesList OBJECT
## ---------------------------------------------------------------------

## A GRangesList object containing exons grouped by transcripts gives
## the same result as above:
exbytx <- exonsBy(txdb, by="tx", use.names=TRUE)
tx_seqs2 <- extractTranscriptsFromGenome(Hsapiens, exbytx)
stopifnot(identical(as.character(tx_seqs2), as.character(tx_seqs1)))

## A sanity check:
stopifnot(identical(unname(sapply(width(exbytx), sum)), width(tx_seqs2)))

## CDSs grouped by transcripts (this extracts only the translated parts
## of the transcripts):
cds_seqs <- extractTranscriptsFromGenome(Hsapiens, cdsBy(txdb, by="tx"))
translate(cds_seqs)

## ---------------------------------------------------------------------
## C. GOING FROM TRANSCRIPT-BASED TO REFERENCE-BASED LOCATIONS
## ---------------------------------------------------------------------
## Get the reference-based locations of the first 4 (5' end)
## and last 4 (3' end) nucleotides in each transcript:
tlocs <- lapply(width(tx_seqs2), function(w) c(1:4, (w-3):w))
tx_strand <- sapply(strand(exbytx), runValue)
## Note that, because of how we made them, 'tlocs', 'start(exbytx)',
## 'end(exbytx)' and 'tx_strand' have the same length, and, for any
## valid positional index, elements at this position are corresponding
## to each other. This is how transcriptLocs2refLocs() expects them
## to be!
rlocs <- transcriptLocs2refLocs(tlocs, start(exbytx), end(exbytx),
             tx_strand, decreasing.rank.on.minus.strand=TRUE)

## ---------------------------------------------------------------------
## D. EXTRACTING WORM TRANSCRIPTS ZC101.3 AND F37B1.1
## ---------------------------------------------------------------------

## Transcript ZC101.3 (is on + strand):
##   Exons starts/ends relative to transcript:
rstarts1 <- c(1, 488, 654, 996, 1365, 1712, 2163, 2453)
rends1 <- c(137, 578, 889, 1277, 1662, 1870, 2410, 2561)
##   Exons starts/ends relative to chromosome:
starts1 <- 14678410 + rstarts1
ends1 <- 14678410 + rends1

## Transcript F37B1.1 (is on - strand):
##   Exons starts/ends relative to transcript:
rstarts2 <- c(1, 325)
rends2 <- c(139, 815)
##   Exons starts/ends relative to chromosome:
starts2 <- 13611188 - rends2
ends2 <- 13611188 - rstarts2

exon_starts <- list(as.integer(starts1), as.integer(starts2))
exon_ends <- list(as.integer(ends1), as.integer(ends2))

library(BSgenome.Celegans.UCSC.ce2)
## Both transcripts are on chrII:
chrII <- Celegans$chrII
tx_seqs <- extractTranscripts(chrII,
                              exonStarts=exon_starts,
                              exonEnds=exon_ends,
                              strand=c("+","-"))

## Same as 'width(tx_seqs)':
transcriptWidths(exonStarts=exon_starts, exonEnds=exon_ends)

transcriptLocs2refLocs(list(c(1:6, 135:140, 1555:1560),
                            c(1:6, 137:142, 625:630)),
                       exonStarts=exon_starts,
                       exonEnds=exon_ends,
                       strand=c("+","-"))

## A sanity check:
ref_locs <- transcriptLocs2refLocs(list(1:1560, 1:630),
                                   exonStarts=exon_starts,
                                   exonEnds=exon_ends,
                                   strand=c("+","-"))
stopifnot(chrII[ref_locs[[1]]] == tx_seqs[[1]])
stopifnot(complement(chrII)[ref_locs[[2]]] == tx_seqs[[2]])

## ---------------------------------------------------------------------
## E. sortExonsByRank()
## ---------------------------------------------------------------------
## Typically used to reorder by decreasing rank the exons in transcripts
## located on the minus strand:
exbytx3 <- sortExonsByRank(exbytx, decreasing.rank.on.minus.strand=TRUE)
exbytx3
tx_seqs3 <- extractTranscriptsFromGenome(Hsapiens, exbytx3,
                            decreasing.rank.on.minus.strand=TRUE)
stopifnot(identical(as.character(tx_seqs3), as.character(tx_seqs1)))
}

\keyword{manip}
