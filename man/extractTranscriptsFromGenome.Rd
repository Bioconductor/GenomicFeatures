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

\author{
  H. Pages
}

\seealso{
  \code{\link[BSgenome:available.genomes]{available.genomes}},
  \code{\link[GenomicFeatures.Hsapiens.UCSC.hg18]{geneHuman}}
}

\examples{
  if (interactive()) {
    library(BSgenome.Hsapiens.UCSC.hg18)  # load the genome
    library(GenomicFeatures.Hsapiens.UCSC.hg18)  # load the gene table
    ## Takes about 30 sec.
    transcripts <- extractTranscriptsFromGenome(Hsapiens, geneHuman())
  }
}
