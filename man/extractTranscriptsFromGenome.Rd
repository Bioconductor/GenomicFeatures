\name{extractTranscriptsFromGenome}

\alias{extractTranscriptsFromGenome}

\title{Extract transcripts from a genome}

\description{
  \code{extractTranscriptsFromGenome} extracts the transcript
  sequences from a BSgenome data package using transcript
  information (exon boundaries) stored in a "gene table".
}

\usage{
  extractTranscriptsFromGenome(genome, genes)
}

\arguments{
  \item{genome}{
    A \code{\link[BSgenome:BSgenome-class]{BSgenome}} object.
    See the \code{\link[BSgenome]{available.genomes}} function
    in the BSgenome package for how to install a genome.
  }
  \item{genes}{
    A data.frame like that returned by
    \code{\link[GenomicFeatures.Hsapiens.UCSC.hg18]{geneHuman}}.
  }
}

\value{
  A \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} object.
}

\note{
  \code{extractTranscriptsFromGenome} is based on the
  \code{\link[Biostrings:extractTranscripts]{extractTranscripts}}
  function defined in the Biostrings package.
  See \code{?`\link[Biostrings:extractTranscripts]{extractTranscripts}`}
  for more information and related functions like
  \code{\link[Biostrings:extractTranscripts]{transcriptLocs2refLocs}}
  for converting transcript-based locations into chromosome-based
  (aka reference-based) locations.
}

\author{
  H. Pages
}

\seealso{
  \code{\link[BSgenome:available.genomes]{available.genomes}},
  \code{\link[GenomicFeatures.Hsapiens.UCSC.hg18]{geneHuman}},
  \code{\link[Biostrings:extractTranscripts]{transcriptLocs2refLocs}}
}

\examples{
  library(GenomicFeatures.Hsapiens.UCSC.hg18)  # load the gene table
  genes <- geneHuman()
  tw <- transcriptWidths(genes$exonStarts, genes$exonEnds)

  if (interactive()) {
    library(BSgenome.Hsapiens.UCSC.hg18)  # load the genome
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
