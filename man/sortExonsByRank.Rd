\name{sortExonsByRank}

% Old stuff
\alias{sortExonsByRank}

\title{Sort exons by rank}

\description{
  WARNING: Starting with BioC 3.0, \code{sortExonsByRank} is defunct
  with no replacement.

  \code{sortExonsByRank} orders (or reorders) by rank the exons
  stored in a \link[GenomicRanges]{GRangesList} object containing
  exons grouped by transcript.
}

\usage{
# DEFUNCT. No replacement.
sortExonsByRank(x, decreasing.rank.on.minus.strand=FALSE)
}

\arguments{
  \item{x}{
    A \link[GenomicRanges]{GRangesList} object, typically coming from
    \code{exonsBy(... , by="tx")}.
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
}

\value{
  A \link[GenomicRanges]{GRangesList} object with one top-level element per
  transcript. More precisely, the returned object has the same "shape" (i.e.
  same length and same number of elements per top-level element) as the input
  \link[GenomicRanges]{GRangesList} object \code{x}.
}

\author{
  H. Pages
}

\keyword{manip}
