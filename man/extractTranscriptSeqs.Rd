\name{extractTranscriptSeqs}

\alias{extractTranscriptSeqs}
\alias{extractTranscriptSeqs,DNAString-method}
\alias{extractTranscriptSeqs,BSgenome-method}

\title{Extract transcript sequences from chromosomes}

\description{
  \code{extractTranscriptSeqs} is a generic function for extracting
  transcript sequences from an object representing a single chromosome
  (e.g. a \link[Biostrings]{DNAString} object) or a collection of
  chromosomes (e.g. a \link[BSgenome]{BSgenome} object).
}

\usage{
extractTranscriptSeqs(x, transcripts, ...)

\S4method{extractTranscriptSeqs}{DNAString}(x, transcripts, strand="+")

\S4method{extractTranscriptSeqs}{BSgenome}(x, transcripts)
}

\arguments{
  \item{x}{
    A \link[Biostrings]{DNAString} or \link[BSgenome]{BSgenome} object.
  }
  \item{transcripts}{
    An object representing the exon ranges of each transcript to extract.
 
    It must be an \link[IRanges]{RangesList} object when \code{x} is a
    \link[Biostrings]{DNAString} object.

    It must be a \link[GenomicRanges]{GRangesList} or \link{TranscriptDb}
    object when \code{x} is a \link[BSgenome]{BSgenome} object.
    If the latter, it's first turned into a \link[GenomicRanges]{GRangesList}
    object with \code{\link{exonsBy}(transcripts, by="tx", use.names=TRUE)}.
  }
  \item{strand}{
    Only supported when \code{x} is a \link[Biostrings]{DNAString} object.

    Can be an atomic vector, a factor, or an Rle object, in which case it
    indicates the strand of each transcript (i.e. all the exons in a
    transcript are considered to be on the same strand).
    More precisely: it's turned into a factor (or factor-\link[IRanges]{Rle})
    that has the "standard strand levels" (this is done by calling the
    \code{\link[BiocGenerics]{strand}} function on it). Then it's recycled
    to the length of \link[IRanges]{RangesList} object \code{transcripts}
    if needed. In the resulting object, the i-th element is interpreted
    as the strand of all the exons in the i-th transcript.

    \code{strand} can also be a list-like object, in which case it indicates
    the strand of each exon, individually. Thus it must have the same
    \emph{shape} as \link[IRanges]{RangesList} object \code{transcripts}
    (i.e. same length plus \code{strand[[i]]} must have the same length
    as \code{transcripts[[i]]} for all \code{i}).

    \code{strand} can only contain \code{"+"} and/or \code{"-"} values.
    \code{"*"} is not allowed.
  }
}

\value{
  A \link[Biostrings]{DNAStringSet} object \emph{parallel} to
  \code{transcripts}, that is, the i-th element in the returned object
  is the sequence of the i-th transcript in \code{transcripts}.
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item The \code{\link[BSgenome]{available.genomes}} function in the
          \pkg{BSgenome} package for checking avaibility of BSgenome
          data packages (and installing the desired one).
          
    \item The \link[GenomicRanges]{GRangesList} class defined and documented
          in the \pkg{GenomicRanges} package.

    \item The \link[IRanges]{RangesList} class defined and documented
          in the \pkg{IRanges} package.

    \item The \code{\link{exonsBy}} function for extracting exon ranges
          grouped by transcript.

    \item The \link[Biostrings]{DNAString} and \link[Biostrings]{DNAStringSet}
          classes defined and documented in the \pkg{Biostrings} package.

    \item The \code{\link[Biostrings]{translate}} function in the
          \pkg{Biostrings} package for translating DNA or RNA sequences
          into amino acid sequences.
  }
}

\examples{
## ---------------------------------------------------------------------
## A FIRST EXAMPLE
## ---------------------------------------------------------------------

## Load a genome:
library(BSgenome.Hsapiens.UCSC.hg18)
genome <- BSgenome.Hsapiens.UCSC.hg18

## Load a TranscriptDb object:
txdb_file <- system.file("extdata", "UCSC_knownGene_sample.sqlite",
                         package="GenomicFeatures")
txdb <- loadDb(txdb_file)

## Check that 'txdb' is based on the hg18 assembly:
txdb

## Extract the exon ranges grouped by transcript from 'txdb':
transcripts <- exonsBy(txdb, by="tx", use.names=TRUE)

## Extract the transcript sequences from the genome:
tx_seqs <- extractTranscriptSeqs(genome, transcripts)
tx_seqs

## A sanity check:
stopifnot(identical(width(tx_seqs), unname(sum(width(transcripts)))))

## ---------------------------------------------------------------------
## USING extractTranscriptSeqs() TO EXTRACT CDS SEQUENCES
## ---------------------------------------------------------------------

cds <- cdsBy(txdb, by="tx", use.names=TRUE)
cds_seqs <- extractTranscriptSeqs(genome, cds)
cds_seqs

## A sanity check:
stopifnot(identical(width(cds_seqs), unname(sum(width(cds)))))

## Note that, alternatively, the CDS sequences can be obtained from the
## transcript sequences by removing the 5' and 3' UTRs:
five_utr_width <- sum(width(fiveUTRsByTranscript(txdb, use.names=TRUE)))
five_utr_width <- five_utr_width[names(cds_seqs)]
five_utr_width[is.na(five_utr_width)] <- 0L
three_utr_width <- sum(width(threeUTRsByTranscript(txdb, use.names=TRUE)))
three_utr_width <- three_utr_width[names(cds_seqs)]
three_utr_width[is.na(three_utr_width)] <- 0L
cds_seqs2 <- narrow(tx_seqs[names(cds_seqs)],
                    start=five_utr_width+1L,
                    end=-(three_utr_width+1L))
stopifnot(identical(as.character(cds_seqs), as.character(cds_seqs2)))

## ---------------------------------------------------------------------
## TRANSLATE THE CDS SEQUENCES
## ---------------------------------------------------------------------

prot_seqs <- translate(cds_seqs, if.fuzzy.codon="solve")

## Note that, by default, translate() uses The Standard Genetic Code to
## translate codons into amino acids. However, depending on the organism,
## a different genetic code might be needed to translate CDS sequences
## located on the mitochodrial chromosome. For example, for vertebrates,
## the following code could be used to correct 'prot_seqs':
SGC1 <- getGeneticCode("SGC1")
chrM_idx <- which(all(seqnames(cds) == "chrM"))
prot_seqs[chrM_idx] <- translate(cds_seqs[chrM_idx], genetic.code=SGC1,
                                 if.fuzzy.codon="solve")
}

\keyword{manip}
