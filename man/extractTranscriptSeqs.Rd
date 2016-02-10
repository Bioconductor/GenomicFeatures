\name{extractTranscriptSeqs}

\alias{extractTranscriptSeqs}
\alias{extractTranscriptSeqs,DNAString-method}
\alias{extractTranscriptSeqs,ANY-method}


\title{Extract transcript (or CDS) sequences from chromosome sequences}

\description{
  \code{extractTranscriptSeqs} extracts transcript (or CDS) sequences from
  an object representing a single chromosome or a collection of chromosomes.
}

\usage{
extractTranscriptSeqs(x, transcripts, ...)

\S4method{extractTranscriptSeqs}{DNAString}(x, transcripts, strand="+")

\S4method{extractTranscriptSeqs}{ANY}(x, transcripts, ...)
}

\arguments{
  \item{x}{
    An object representing a single chromosome or a collection of chromosomes.
    More precisely, \code{x} can be a \link[Biostrings]{DNAString} object
    (single chromosome), or a \link[BSgenome]{BSgenome} object (collection
    of chromosomes).

    Other objects representing a collection of chromosomes are supported
    (e.g. \link[Rsamtools]{FaFile} objects in the \pkg{Rsamtools} package)
    as long as \code{\link[GenomeInfoDb]{seqinfo}} and
    \code{\link[Biostrings]{getSeq}} work on them.
  }
  \item{transcripts}{
    An object representing the exon ranges of each transcript to extract.

    More precisely:
    \itemize{
      \item If \code{x} is a \link[Biostrings]{DNAString} object, then
            \code{transcripts} must be an \link[IRanges]{RangesList} object.

      \item If \code{x} is a \link[BSgenome]{BSgenome} object or any object
            representing a collection of chromosomes, then \code{transcripts}
            must be a \link[GenomicRanges]{GRangesList} object or any object
            for which \code{\link{exonsBy}} is implemented (e.g. a \link{TxDb}
            object). If the latter, then it's first turned into a
            \link[GenomicRanges]{GRangesList} object with
            \code{\link{exonsBy}(transcripts, by="tx", ...)}.
    }

    Note that, for each transcript, the exons must be ordered by ascending
    \emph{rank}, that is, by ascending position \emph{in the transcript}
    (when going in the 5' to 3' direction). This generally means (but not
    always) that they are also ordered from 5' to 3' on the reference genome.
    More precisely:
    \itemize{
      \item For a transcript located on the plus strand, the exons will
            typically (but not necessarily) be ordered by ascending position
            on the reference genome.
      \item For a transcript located on the minus strand, the exons will
            typically (but not necessarily) be ordered by descending position
            on the reference genome.
    }
    If \code{transcripts} was obtained with \code{\link{exonsBy}} (see above),
    then the exons are guaranteed to be ordered by ascending rank. See
    \code{?\link{exonsBy}} for more information.
  }
  \item{...}{
    Additional arguments, for use in specific methods.

    For the default method, additional arguments are allowed only when
    \code{transcripts} is not a \link[GenomicRanges]{GRangesList} object,
    in which case they are passed to the internal call to \code{\link{exonsBy}}
    (see above).
  }
  \item{strand}{
    Only supported when \code{x} is a \link[Biostrings]{DNAString} object.

    Can be an atomic vector, a factor, or an \link[S4Vectors]{Rle} object,
    in which case it indicates the strand of each transcript (i.e. all the
    exons in a transcript are considered to be on the same strand).
    More precisely: it's turned into a factor (or factor-\link[S4Vectors]{Rle})
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
  \code{transcripts}, that is, the i-th element in it is the sequence
  of the i-th transcript in \code{transcripts}.
}

\author{
  H. Pages
}

\seealso{
  \itemize{
    \item \code{\link{coverageByTranscript}} for computing coverage by
          transcript (or CDS) of a set of ranges.

    \item \code{\link{transcriptLengths}} for extracting the transcript
          lengths from a \link{TxDb} object.

    \item The \code{\link{transcriptLocs2refLocs}} function for converting
          transcript-based locations into reference-based locations.

    \item The \code{\link[BSgenome]{available.genomes}} function in the
          \pkg{BSgenome} package for checking avaibility of BSgenome
          data packages (and installing the desired one).
          
    \item The \link[Biostrings]{DNAString} and \link[Biostrings]{DNAStringSet}
          classes defined and documented in the \pkg{Biostrings} package.

    \item The \code{\link[Biostrings]{translate}} function in the
          \pkg{Biostrings} package for translating DNA or RNA sequences
          into amino acid sequences.

    \item The \link[GenomicRanges]{GRangesList} class defined and documented
          in the \pkg{GenomicRanges} package.

    \item The \link[IRanges]{RangesList} class defined and documented
          in the \pkg{IRanges} package.

    \item The \code{\link{exonsBy}} function for extracting exon ranges
          grouped by transcript.

    \item The \link{TxDb} class.
  }
}

\examples{
## ---------------------------------------------------------------------
## 1. A TOY EXAMPLE
## ---------------------------------------------------------------------

library(Biostrings)

## A chromosome of length 30:
x <- DNAString("ATTTAGGACACTCCCTGAGGACAAGACCCC")

## 2 transcripts on 'x':
tx1 <- IRanges(1, 8)            # 1 exon
tx2 <- c(tx1, IRanges(12, 30))  # 2 exons
transcripts <- IRangesList(tx1=tx1, tx2=tx2)
extractTranscriptSeqs(x, transcripts)

## By default, all the exons are considered to be on the plus strand.
## We can use the 'strand' argument to tell extractTranscriptSeqs()
## to extract them from the minus strand.

## Extract all the exons from the minus strand:
extractTranscriptSeqs(x, transcripts, strand="-")

## Note that, for a transcript located on the minus strand, the exons
## should typically be ordered by descending position on the reference
## genome in order to reflect their rank in the transcript:
extractTranscriptSeqs(x, IRangesList(tx1=tx1, tx2=rev(tx2)), strand="-")

## Extract the exon of the 1st transcript from the minus strand:
extractTranscriptSeqs(x, transcripts, strand=c("-", "+"))

## Extract the 2nd exon of the 2nd transcript from the minus strand:
extractTranscriptSeqs(x, transcripts, strand=list("-", c("+", "-")))

## ---------------------------------------------------------------------
## 2. A REAL EXAMPLE
## ---------------------------------------------------------------------

## Load a genome:
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19

## Load a TxDb object:
txdb_file <- system.file("extdata", "hg19_knownGene_sample.sqlite",
                         package="GenomicFeatures")
txdb <- loadDb(txdb_file)

## Check that 'txdb' is based on the hg19 assembly:
txdb

## Extract the exon ranges grouped by transcript from 'txdb':
transcripts <- exonsBy(txdb, by="tx", use.names=TRUE)

## Extract the transcript sequences from the genome:
tx_seqs <- extractTranscriptSeqs(genome, transcripts)
tx_seqs

## A sanity check:
stopifnot(identical(width(tx_seqs), unname(sum(width(transcripts)))))

## Note that 'tx_seqs' can also be obtained with:
extractTranscriptSeqs(genome, txdb, use.names=TRUE)

## ---------------------------------------------------------------------
## 3. USING extractTranscriptSeqs() TO EXTRACT CDS SEQUENCES
## ---------------------------------------------------------------------

cds <- cdsBy(txdb, by="tx", use.names=TRUE)
cds_seqs <- extractTranscriptSeqs(genome, cds)
cds_seqs

## A sanity check:
stopifnot(identical(width(cds_seqs), unname(sum(width(cds)))))

## Note that, alternatively, the CDS sequences can be obtained from the
## transcript sequences by removing the 5' and 3' UTRs:
tx_lens <- transcriptLengths(txdb, with.utr5_len=TRUE, with.utr3_len=TRUE)
stopifnot(identical(tx_lens$tx_name, names(tx_seqs)))  # sanity
## Keep the rows in 'tx_lens' that correspond to a sequence in 'cds_seqs'
## and put them in the same order as in 'cds_seqs':
m <- match(names(cds_seqs), names(tx_seqs))
tx_lens <- tx_lens[m, ]
utr5_width <- tx_lens$utr5_len
utr3_width <- tx_lens$utr3_len
cds_seqs2 <- narrow(tx_seqs[m],
                    start=utr5_width+1L, end=-(utr3_width+1L))
stopifnot(identical(as.character(cds_seqs2), as.character(cds_seqs)))

## ---------------------------------------------------------------------
## 4. TRANSLATE THE CDS SEQUENCES
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
