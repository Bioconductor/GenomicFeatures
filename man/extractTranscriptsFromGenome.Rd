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
    object, or a data frame like that returned by
    \code{\link[GenomicFeatures.Hsapiens.UCSC.hg18]{geneHuman}}.
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
  A \link[Biostrings]{DNAStringSet} object.
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
  \code{\link[GenomicFeatures.Hsapiens.UCSC.hg18]{geneHuman}},
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
  transcripts <- extractTranscriptsFromGenome(Hsapiens, txdb)
  transcripts

  ## ---------------------------------------------------------------------
  ## B. USING A GRangesList OBJECT
  ## ---------------------------------------------------------------------
  ## Exons grouped by transcripts (gives the same result as above except
  ## that now transcripts are named by their internal id i.e. by tx_id
  ## instead of tx_name):
  extractTranscriptsFromGenome(Hsapiens, exonsBy(txdb))
  ## CDSs grouped by transcripts (this extracts only the translated parts
  ## of the transcripts):
  cds <- extractTranscriptsFromGenome(Hsapiens, cdsBy(txdb))

  ## ---------------------------------------------------------------------
  ## C. USING A UCSC-LIKE DATA FRAME
  ## ---------------------------------------------------------------------
  ## IMPORTANT NOTE: This is provided for compatibility with the old
  ## GenomicFeatures.*.UCSC.* packages and might be removed at any time.
  library(GenomicFeatures.Hsapiens.UCSC.hg18)  # load the gene table
  genes <- geneHuman()
  library(Biostrings)  # for transcriptWidths()
  tw <- transcriptWidths(genes$exonStarts, genes$exonEnds)

  if (interactive()) {
    ## Takes about 30 sec.:
    transcripts <- extractTranscriptsFromGenome(Hsapiens, genes)
    ## Sanity check:
    stopifnot(identical(width(transcripts), tw))
  }

  ## Get the reference-based locations of the first 4 (5' end)
  ## and last 4 (3' end) nucleotides in each transcript:
  tlocs <- lapply(tw, function(w) c(1:4, (w-3):w))
  rlocs <- transcriptLocs2refLocs(tlocs, genes$exonStarts, genes$exonEnds,
               genes$strand, reorder.exons.on.minus.strand=TRUE)
  
}
