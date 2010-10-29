\name{extractTranscriptsFromGenome}

\alias{extractTranscriptsFromGenome}

\title{Extract transcripts from a genome}

\description{
  \code{extractTranscriptsFromGenome} extracts the transcript
  sequences from a BSgenome data package using transcript
  information (exon boundaries) stored in a "gene table".
}

\usage{
extractTranscriptsFromGenome(genome, txdb, use.names=TRUE)
}

\arguments{
  \item{genome}{
    A \code{\link[BSgenome]{BSgenome}} object.
    See the \code{\link[BSgenome]{available.genomes}} function
    in the BSgenome package for how to install a genome.
  }
  \item{txdb}{
    A \link{TranscriptDb} object, a \link[GenomicRanges]{GRangesList}
    object, or a UCSC-style data frame (i.e. 1 row per transcript and
    it must have at least the \code{"name"}, \code{"chrom"}, \code{"strand"},
    \code{"exonStarts"} and \code{"exonEnds"} columns where each value
    in exonStarts or exonEnds is a comma-separated list of integers --
    note that unlike what UCSC does the start values here must be 1-based,
    not 0-based).
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
}

\value{
  A named \link[Biostrings]{DNAStringSet} object.

  Element positions are preserved with respect to \code{txdb} when
  \code{txdb} is not a \link{TranscriptDb} object.
  More precisely, if \code{txdb} is a \link[GenomicRanges]{GRangesList}
  then elements in the output correspond to elements in \code{txdb}
  and the names are preserved.
  If \code{txdb} is a UCSC-style data frame then elements in the output
  correspond to rows in \code{txdb} and the names are taken from the
  \code{"name"} column.
}

\note{
  \code{extractTranscriptsFromGenome} is based on the
  \code{\link[Biostrings]{extractTranscripts}}
  function defined in the Biostrings package.
  See \code{?`\link[Biostrings]{extractTranscripts}`}
  for more information and related functions like
  \code{\link[Biostrings]{transcriptLocs2refLocs}}
  for converting transcript-based locations into chromosome-based
  (aka reference-based) locations.
}

\author{
  H. Pages
}

\seealso{
  \code{\link[BSgenome]{available.genomes}},
  \code{\link[Biostrings]{transcriptLocs2refLocs}}
}

\examples{
  library(BSgenome.Hsapiens.UCSC.hg18)  # load the genome

  ## ---------------------------------------------------------------------
  ## A. USING A TranscriptDb OBJECT
  ## ---------------------------------------------------------------------
  txdb_file <- system.file("extdata", "UCSC_knownGene_sample.sqlite",
                           package="GenomicFeatures")
  txdb <- loadFeatures(txdb_file)
  txseqs <- extractTranscriptsFromGenome(Hsapiens, txdb)
  txseqs

  ## ---------------------------------------------------------------------
  ## B. USING A GRangesList OBJECT
  ## ---------------------------------------------------------------------
  ## Exons grouped by transcripts (gives the same result as above except
  ## that now transcripts are named by their internal id i.e. by tx_id
  ## instead of tx_name):
  exbytx <- exonsBy(txdb)
  txseqs2 <- extractTranscriptsFromGenome(Hsapiens, exbytx)
  ## A sanity check:
  stopifnot(identical(unname(sapply(width(exbytx), sum)), width(txseqs2)))

  ## CDSs grouped by transcripts (this extracts only the translated parts
  ## of the transcripts):
  cds <- extractTranscriptsFromGenome(Hsapiens, cdsBy(txdb))

  ## ---------------------------------------------------------------------
  ## C. GOING FROM TRANSCRIPT-BASED TO REFERENCE-BASED LOCATIONS
  ## ---------------------------------------------------------------------
  ## Get the reference-based locations of the first 4 (5' end)
  ## and last 4 (3' end) nucleotides in each transcript:
  tlocs <- lapply(width(txseqs2), function(w) c(1:4, (w-3):w))
  tx_strand <- sapply(strand(exbytx), runValue)
  ## Note that, because of how we made them, 'tlocs', 'start(exbytx)',
  ## 'end(exbytx)' and 'tx_strand' have the same length, and, for any
  ## valid positional index, elements at this position are corresponding
  ## to each other. This is how transcriptLocs2refLocs() expects them
  ## to be!
  rlocs <- transcriptLocs2refLocs(tlocs, start(exbytx), end(exbytx),
               tx_strand, reorder.exons.on.minus.strand=TRUE)
}

\keyword{manip}
