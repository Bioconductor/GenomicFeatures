\name{mapToTranscripts}
\alias{coordinate-mapping}

\alias{mapToTranscripts}
\alias{mapToTranscripts,GenomicRanges,GenomicRanges-method}
\alias{mapToTranscripts,GenomicRanges,GRangesList-method}
\alias{mapToTranscripts,ANY,TxDb-method}
\alias{pmapToTranscripts}
\alias{pmapToTranscripts,GenomicRanges,GenomicRanges-method}
\alias{pmapToTranscripts,GenomicRanges,GRangesList-method}
\alias{pmapToTranscripts,GRangesList,GRangesList-method}

\alias{mapFromTranscripts}
\alias{mapFromTranscripts,GenomicRanges,GenomicRanges-method}
\alias{mapFromTranscripts,GenomicRanges,GRangesList-method}
\alias{pmapFromTranscripts}
\alias{pmapFromTranscripts,Ranges,GenomicRanges-method}
\alias{pmapFromTranscripts,Ranges,GRangesList-method}
\alias{pmapFromTranscripts,GenomicRanges,GenomicRanges-method}
\alias{pmapFromTranscripts,GenomicRanges,GRangesList-method}

\title{Map range coordinates between transcripts and genome space}

\description{
  Map range coordinates between features in the transcriptome and
  genome (reference) space.

  See \code{?\link[GenomicAlignments]{mapToAlignments}} in the 
  \pkg{GenomicAlignments} package for mapping coordinates between 
  reads (local) and genome (reference) space using a CIGAR alignment.
}

\usage{
## mapping to transcripts
\S4method{mapToTranscripts}{GenomicRanges,GenomicRanges}(x, transcripts, 
          ignore.strand = FALSE) 
\S4method{mapToTranscripts}{GenomicRanges,GRangesList}(x, transcripts, 
          ignore.strand = FALSE, intronJunctions=FALSE) 
\S4method{mapToTranscripts}{ANY,TxDb}(x, transcripts, ignore.strand = FALSE, 
          extractor.fun = GenomicFeatures::transcripts, ...) 
\S4method{pmapToTranscripts}{GenomicRanges,GRangesList}(x, transcripts, 
          ignore.strand = FALSE) 

## mapping from transcripts
\S4method{mapFromTranscripts}{GenomicRanges,GRangesList}(x, transcripts, 
          ignore.strand = FALSE) 
\S4method{pmapFromTranscripts}{GenomicRanges,GRangesList}(x, transcripts, 
          ignore.strand = FALSE)
\S4method{pmapFromTranscripts}{Ranges,GRangesList}(x, transcripts)
}

\arguments{
  \item{x}{
    \code{\linkS4class{GenomicRanges}} object of positions to be mapped.
    The seqnames of \code{x} are used in \code{mapFromTranscripts}, i.e.,
    when mapping from transcripts to the genome. In the case of
    \code{pmapFromTranscripts}, \code{x} can be a \code{Ranges} object.
  }
  \item{transcripts}{
    A named \code{\linkS4class{GenomicRanges}} or
    \code{\linkS4class{GRangesList}} object used to map between \code{x} and the
    result. The ranges can be any feature in the transcriptome extracted from a
    \code{TxDb} (e.g., introns, exons, cds regions). See ?\code{transcripts} and
    ?\code{transcriptsBy} for a list of extractor functions. 
 
    The \code{transcripts} object must have names. When mapping from
    transcripts to the genome, they are used to determine mapping pairs; in
    the reverse direction they become the seqlevels of the output object.
  }
  \item{ignore.strand}{
    When \code{ignore.strand} is TRUE, strand is ignored in overlaps operations 
    (i.e., all strands are considered "+") and the strand in the output is '*'.

    When \code{ignore.strand} is FALSE strand in the output is taken from the 
    \code{transcripts} argument. When \code{transcripts} is a
    \code{GRangesList}, all inner list elements of a common list element must 
    have the same strand or an error is thrown.

    Mapped position is computed by counting from the transcription start site 
    (TSS) and is not affected by the value of \code{ignore.strand}.
  }
  \item{intronJunctions}{
    Logical to indicate if intron ranges in \code{x} should return a range that
    spans the intron gap. Ranges in \code{x} that do not fall \sQuote{within} a
    range in \code{transcripts} are considered a non-hit and have no mapping.
    When \code{intronJunctions} is TRUE, these no-hits map to a range of width
    2 that spans the junction (i.e., the start and end of the ranges they fall
    between). The range in \code{x} must fall between ranges in
    \code{transcripts}. Ranges that 'fall off' either end or overlap both the
    gap and a range in \code{transcripts} are not considered intron ranges and
    will not be mapped.

    This argument is only supported in \code{mapToTranscripts} when 
    \code{transcripts} is a GRangesList.
  }
  \item{extractor.fun}{
    Function to extract genomic features from a \code{TxDb} object.

    This argument is only applicable to \code{mapToTranscripts} when 
    \code{transcripts} is a \code{TxDb} object. The \code{extractor} should
    be the name of a function (not a character()) described on the 
    \code{?transcripts}, \code{?transcriptsBy}, or \code{?microRNAs} man page.

    Valid \code{extractor} functions:
    \itemize{
      \item{transcripts}  ## default
      \item{exons}
      \item{cds}
      \item{genes}
      \item{promoters}
      \item{disjointExons}
      \item{transcriptsBy}
      \item{exonsBy}
      \item{cdsBy}
      \item{intronsByTranscript}
      \item{fiveUTRsByTranscript}
      \item{threeUTRsByTranscript}
      \item{microRNAs}
      \item{tRNAs}
    }    
  }
  \item{\dots}{
    Additional arguments passed to \code{extractor.fun} functions. 
  }
}

\details{

  In GenomicFeatures >= 1.21.10, the default for \code{ignore.strand} was
  changed to \code{FALSE} for consistency with other methods in GenomicRanges
  and GenomicAlignments. Additionally, the mapped position is computed
  from the TSS and does not depend on the \code{ignore.strand} argument.
  See the section on \code{ignore.strand} for details.
  
  \itemize{
    \item{\code{mapToTranscripts}, \code{pmapToTranscripts}}{
      The genomic range in \code{x} is mapped to the local position in the
      \code{transcripts} ranges. A successful mapping occurs when \code{x}
      is completely within the \code{transcripts} range, equivalent to: 
      \preformatted{
      findOverlaps(..., type="within")
      }
      Transcriptome-based coordinates start counting at 1 at the beginning
      of the \code{transcripts} range and return positions where \code{x}
      was aligned. The seqlevels of the return object are taken from the
      \code{transcripts} object and should be transcript names. In this
      direction, mapping is attempted between all elements of \code{x} and
      all elements of \code{transcripts}.

      \code{mapToTranscripts} uses \code{findOverlaps} to map ranges in
      \code{x} to ranges in \code{transcripts}. This method does not return 
      unmapped ranges.

      \code{pmapToTranscripts} maps the i-th range in \code{x} to the
      i-th range in \code{transcripts}. Recycling is supported for both
      \code{x} and \code{transcripts} when either is length == 1L; otherwise
      the lengths must match. Ranges in \code{x} that do not map (out of bounds 
      or strand mismatch) are returned as zero-width ranges starting at 0. 
      These ranges are given the seqname of "UNMAPPED". 
    }
    \item{\code{mapFromTranscripts}, \code{pmapFromTranscripts}}{
      The transcript-based position in \code{x} is mapped to genomic coordinates
      using the ranges in \code{transcripts}. A successful mapping occurs when
      the following is TRUE:
      \preformatted{
      width(transcripts) >= start(x) + width(x)
      }
      \code{x} is aligned to \code{transcripts} by moving in \code{start(x)}
      positions in from the beginning of the \code{transcripts} range.  The
      seqlevels of the return object are chromosome names.

      \code{mapFromTranscripts} uses the seqname of \code{x} and the names 
      of \code{transcripts} to determine mapping pairs (vs attempting to match 
      all possible pairs). Name matching is motivated by use cases such as 
      differentially expressed regions where the expressed regions in \code{x} 
      would only be related to a subset of regions in \code{transcripts}. 
      This method does not return unmapped ranges.
     
      \code{pmapFromTranscripts} maps the i-th range in \code{x} to the i-th
      range in \code{transcripts} and therefore does not use name matching.
      Recycling is supported in \code{pmapFromTranscripts} when either
      \code{x} or \code{transcripts} is length == 1L; otherwise the lengths 
      must match. Ranges in \code{x} that do not map (out of bounds or strand 
      mismatch) are returned as zero-width ranges starting at 0. These ranges 
      are given the seqname of "UNMAPPED". 
    }
  }
}

\value{

  \code{pmapToTranscripts} returns a \code{GRanges} the same length as
  \code{x}.
  
  \code{pmapFromTranscripts} returns a \code{GRanges} when \code{transcripts}
  is a \code{GRanges} and a \code{GRangesList} when \code{transcripts}
  is a \code{GRangesList}. In both cases the return object is the same
  length as \code{x}. The rational for returning the \code{GRangesList} is 
  to preserve exon structure; ranges in a list element that are not overlapped 
  by \code{x} are returned as a zero-width range. The \code{GRangesList}
  return object will have no seqlevels called "UNMAPPED"; those will only
  occur when a \code{GRanges} is returned.

  \code{mapToTranscripts} and \code{mapFromTranscripts} return \code{GRanges}
  objects that vary in length similar to a \code{Hits} object. The result 
  contains mapped records only; strand mismatch and out of bound ranges are 
  not returned. \code{xHits} and \code{transcriptsHits} metadata columns 
  (similar to the \code{queryHits} and \code{subjectHits} of a \code{Hits}
  object) indicate elements of \code{x} and \code{transcripts} used in 
  the mapping.

  When mapping to transcript coordinates, seqlevels of the output are the names 
  on the \code{transcripts} object and most often these will be transcript
  names. When mapping to the genome, seqlevels of the output are the seqlevels 
  of \code{transcripts} which are usually chromosome names.
}

\seealso{
  \itemize{
    \item \code{?\link[GenomicAlignments]{mapToAlignments}} in the 
    \pkg{GenomicAlignments} package for methods mapping between
    reads and genome space using a CIGAR alignment.
  }
}

\author{V. Obenchain, M. Lawrence and H. Pagès}

\examples{

## ---------------------------------------------------------------------
## A. Basic Use
## ---------------------------------------------------------------------

## ------------------------------------
## (i) Map from genome to transcript:

## The seqnames of the output are the transcript names, not chromosomes. For
## this reason 'transcripts' must be named. 
x <- GRanges("A", IRanges(16, 18))
gr1 <- GRanges("A", IRanges(1, 10, names="tx_a"))
gr2 <- GRanges("A", IRanges(15, 20, names="tx_b"))

## 'transcripts' as GRanges:
mapToTranscripts(x, gr2)

## 'transcripts' as GRangesList:
mapToTranscripts(x, GRangesList("tx_c" = c(gr1, gr2)))

## Round trip from genomic -> transcript -> genomic coordinates:
tx_coord <- mapToTranscripts(x, gr2)
mapFromTranscripts(tx_coord, gr2)

## ------------------------------------
## (ii) Map from transcript to genome:

## A prerequisite for mapping from transcript -> genome is that the seqname 
## of the range in 'x' match the name of the range in 'transcripts'. Here 
## the seqname of 'x' is "TX_1" and mapping is only attempted with the second 
## range in 'gr':
x <- GRanges("TX_1", IRanges(5, 10))
gr <- GRanges("chr3", IRanges(c(1, 1), width=50, names=c("TX_2", "TX_1")))
mapFromTranscripts(x, gr)

## ------------------------------------
## (iii) Element-wise versions:

## Recycling is supported when length(transcripts) == 1; otherwise the
## lengths of 'x' and 'transcripts' must be the same.
x <- GRanges("A", IRanges(c(1, 5, 10), width=1))
transcripts <- GRanges("A", IRanges(4, 7))
pmapToTranscripts(x, transcripts)


## ---------------------------------------------------------------------
## B. Map local sequence locations to the genome
## ---------------------------------------------------------------------

## NAGNAG alternative splicing plays an essential role in biological processes 
## and represents a highly adaptable system for posttranslational regulation 
## of gene function. The majority of NAGNAG studies largely focus on messenger 
## RNA. A study by Sun, Lin, and Yan 
## (http://www.hindawi.com/journals/bmri/2014/736798/) demonstrated that
## NAGNAG splicing is also operative in large intergenic noncoding RNA
## (lincRNA). 

## One finding of interest was that linc-POLR3G-10 exhibited two NAGNAG 
## acceptors located in two distinct transcripts: TCONS_00010012 and 
## TCONS_00010010. 

## Extract the exon coordinates of TCONS_00010012 and TCONS_00010010: 
lincrna <- c("TCONS_00010012", "TCONS_00010010")
library(TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts)
txdb <- TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts
exons <- exonsBy(txdb, by="tx", use.names=TRUE)[lincrna]
exons

## The two NAGNAG acceptors were identified in the upstream region of 
## the fourth and fifth exons located in TCONS_00010012.
## Extract the sequences for transcript TCONS_00010012:
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19
exons_seq <- getSeq(genome, exons[[1]])

## TCONS_00010012 has 4 exons:
exons_seq

## The most common triplet among the lincRNA sequences was CAG. Identify
## the location of this pattern in all exons.
cag_loc <- vmatchPattern("CAG", exons_seq)

## Convert the first occurance of CAG in each exon back to genome coordinates.
first_loc <- do.call(c, sapply(cag_loc, "[", 1, simplify=TRUE))
pmapFromTranscripts(first_loc, exons[[1]])


## -----------------------------------------------------------------------
## C. Map 3'UTR variants to genome coordinates
## -----------------------------------------------------------------------

## A study by Skeeles et. al (PLoS ONE 8(3): e58609. doi:
## 10.1371/journal.pone.0058609) investigated the impact of 3'UTR variants 
## on the expression of cancer susceptibility genes.

## 8 candidate miRNA genes on chromosome 12 were used to test for 
## differential luciferase expression in mice. In Table 2 of the manuscript
## variant locations are given as nucleotide position within the gene.
geneNames <- c("Bcap29", "Dgkb", "Etv1", "Hbp1", "Hbp1", "Ifrd1", 
               "Ifrd1", "Pik3cg", "Pik3cg", "Tspan13", "Twistnb")
starts <- c(1409, 3170, 3132, 2437, 2626, 3239, 3261, 4947, 4979, 958, 1489) 
snps <- GRanges(geneNames, IRanges(starts, width=1))

## To map these transcript-space coordinates to the genome we need gene ranges
## in genome space.
library(org.Mm.eg.db)
geneid <- select(org.Mm.eg.db, unique(geneNames), "ENTREZID", "SYMBOL")
geneid

## Extract the gene regions:
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)[geneid$ENTREZID]

## A prerequesite of the mapping from transcript space to genome space
## is that seqnames in 'x' match names in 'transcripts'. Rename
## 'genes' with the appropriate gene symbol.
names(genes) <- geneid$SYMBOL 

## The xHits and transcriptsHits metadta columns indicate which ranges in
## 'snps' and 'genes' were involved in the mapping.
mapFromTranscripts(snps, genes)


## -----------------------------------------------------------------------
## D. Map dbSNP variants to cds or cDNA coordinates
## -----------------------------------------------------------------------

## The GIPR gene encodes a G-protein coupled receptor for gastric inhibitory 
## polypeptide (GIP). Originally GIP was identified to inhibited gastric acid 
## secretion and gastrin release but was later demonstrated to stimulate 
## insulin release in the presence of elevated glucose.

## In this example 5 SNPs located in the GIPR gene are mapped to cDNA 
## coordinates. A list of SNPs in GIPR can be downloaded from dbSNP or NCBI.
rsids <- c("rs4803846", "rs139322374", "rs7250736", "rs7250754", "rs9749185")

## Extract genomic coordinates with a SNPlocs package.
library(SNPlocs.Hsapiens.dbSNP141.GRCh38)
snps <- snpid2grange(SNPlocs.Hsapiens.dbSNP141.GRCh38, rsids)

## Gene regions of GIPR can be extracted from a TxDb package of compatible
## build. The TxDb package uses Entrez gene identifiers and GIPR is a gene 
## symbol. Conversion between gene symbols and Entrez gene IDs is done by 
## calling select() on an organism db package.
library(org.Hs.eg.db)
geneid <- select(org.Hs.eg.db, "GIPR", "ENTREZID", "SYMBOL")

## The transcriptsBy() extractor returns a range for each transcript that
## includes the UTR and exon regions (i.e., cDNA).
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txbygene <- transcriptsBy(txdb, "gene")
cDNA <- txbygene[geneid$ENTREZID]
cDNA

## Before mapping, the chromosome names (seqlevels) in the two objects must 
## be harmonized. The style for 'snps' is dbSNP and 'cDNA' is UCSC.
seqlevelsStyle(snps)
seqlevelsStyle(cDNA)

## Modify the style and genome in 'snps' to match 'cDNA'.
seqlevelsStyle(snps) <- seqlevelsStyle(cDNA)
genome(snps) <- genome(cDNA)

## The 'cDNA' object is a GRangesList of length 1. This single list element
## contains the cDNA range for 4 different transcripts. To map to each
## transcript individually 'cDNA' must be unlisted before mapping.

## Map all 5 SNPS to all 4 transcripts:
mapToTranscripts(snps, unlist(cDNA))

## Map the first SNP to transcript uc002pct.1 and the second to uc002pcu.1.
pmapToTranscripts(snps[1:2], unlist(cDNA)[1:2])

## The cdsBy() extractor returns coding regions by gene or by transcript.
## Extract the coding regions for transcript uc002pct.1.
cds <- cdsBy(txdb, "tx", use.names=TRUE)["uc002pct.1"]
cds

## The 'cds' object is a GRangesList of length 1 containing all cds ranges
## for the single transcript uc002pct.1.

## To map to the concatenated group of ranges leave 'cds' as a GRangesList.
mapToTranscripts(snps, cds)

## Only the second SNP could be mapped. Unlisting the 'cds' object maps the
## SNPs to the individual cds ranges (vs the concatenated range).
mapToTranscripts(snps[2], unlist(cds))

## The location is the same because the SNP hit the first cds range. If the
## transcript were on the "-" strand the difference in concatenated vs
## non-concatenated position would be more obvious.

## Change the strand of 'cds':
strand(cds) <- "-"

## Re-map using 'ignore.strand'. The 'ignore.strand' argument is used
## in overlaps operations but does not affect the mapped position.

## Map to concatenated group of cds regions:
mapToTranscripts(snps[2], cds, ignore.strand=TRUE)

## Map to individual cds regions:
mapToTranscripts(snps[2], unlist(cds), ignore.strand=TRUE)
}

\keyword{methods}
\keyword{utilities}
