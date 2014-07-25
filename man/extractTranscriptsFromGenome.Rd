\name{extractTranscriptsFromGenome}

% Old stuff
\alias{extractTranscripts}
\alias{extractTranscriptsFromGenome}
\alias{sortExonsByRank}

\title{Various DEFUNCT tools for extracting transcript sequences}

\description{
  WARNING: Starting with BioC 3.0, the tools described in this man page
  are DEFUNCT.

  \code{extractTranscriptsFromGenome} extracts the transcript
  sequences from a BSgenome data package using the transcript
  information (exon boundaries) stored in a \link{TxDb}
  or \link[GenomicRanges]{GRangesList} object.
  DEFUNCT. Please use \code{\link{extractTranscriptSeqs}} instead.

  \code{extractTranscripts} extracts a set of transcripts
  from a single DNA sequence.
  DEFUNCT. Please use \code{\link{extractTranscriptSeqs}} instead.

  Related utility:

  \code{sortExonsByRank} orders (or reorders) by rank the exons
  stored in a \link[GenomicRanges]{GRangesList} object containing
  exons grouped by transcript.
  DEPRECATED. No replacement.
}

\usage{
# DEFUNCT. Please use extractTranscriptSeqs() instead.
extractTranscriptsFromGenome(genome, txdb,
                             decreasing.rank.on.minus.strand=FALSE,
                             use.names=TRUE)

# DEFUNCT. Please use extractTranscriptSeqs() instead.
extractTranscripts(x,
        exonStarts=list(), exonEnds=list(), strand=character(0),
        decreasing.rank.on.minus.strand=FALSE)

# DEPRECATED. No replacement.
sortExonsByRank(x, decreasing.rank.on.minus.strand=FALSE)
}

\arguments{
  \item{genome}{
    A \link[BSgenome]{BSgenome} object.
    See the \code{\link[BSgenome]{available.genomes}} function
    in the BSgenome package for how to install a genome.
  }
  \item{txdb}{
    A \link{TxDb} object or a \link[GenomicRanges]{GRangesList}
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
    a \link{TxDb} object.
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
}

\value{
  For \code{extractTranscriptsFromGenome}: A named
  \link[Biostrings]{DNAStringSet} object with one element per transcript.
  When \code{txdb} is a \link[GenomicRanges]{GRangesList} object, elements
  in the output align with elements in the input (\code{txdb}), and they
  have the same names.

  For \code{extractTranscripts}: A \link[Biostrings]{DNAStringSet} object
  with one element per transcript.

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
  \code{\link{extractTranscriptSeqs}} as a replacement for
  \code{extractTranscriptsFromGenome} and \code{extractTranscripts}.
}

\keyword{manip}
