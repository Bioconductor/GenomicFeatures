\name{mapToTranscripts}
\alias{coordinate-mapping}

\alias{mapToTranscripts}
\alias{mapToTranscripts,GenomicRanges,GenomicRanges-method}
\alias{mapToTranscripts,GenomicRanges,GRangesList-method}

\alias{pmapToTranscripts}
\alias{pmapToTranscripts,Ranges,GenomicRanges-method}
\alias{pmapToTranscripts,GenomicRanges,GenomicRanges-method}
\alias{pmapToTranscripts,GenomicRanges,GRangesList-method}


\title{Map range coordinates between transcripts and genome space}

\description{
  Map range coordinates between features in the transcriptome and
  genome (reference) space.

  See \code{?\link[GenomicAlignments]{mapToAlignments}} in the 
  \pkg{GenomicAlignments} package for mapping coordinates between 
  reads (local) and genome (reference) space using a CIGAR alignment.
}

\usage{
\S4method{mapToTranscripts}{GenomicRanges,GenomicRanges}(x, transcripts, 
          reverse = FALSE, ignore.strand = TRUE, ...) 
\S4method{mapToTranscripts}{GenomicRanges,GRangesList}(x, transcripts, 
          reverse = FALSE, ignore.strand = TRUE, ...) 

\S4method{pmapToTranscripts}{Ranges,GenomicRanges}(x, transcripts, 
          reverse = FALSE, ...) 
\S4method{pmapToTranscripts}{GenomicRanges,GenomicRanges}(x, transcripts, 
          reverse = FALSE, ignore.strand = TRUE, ...) 
\S4method{pmapToTranscripts}{GenomicRanges,GRangesList}(x, transcripts, 
          reverse = FALSE, ignore.strand = TRUE, ...) 
}

\arguments{
  \item{x}{
    \code{\linkS4class{GenomicRanges}} object of positions to be mapped.
    When \code{reverse = TRUE} (i.e., the mapping is to genome space)
    \code{x} must have names.
  }
  \item{transcripts}{
    A named \code{\linkS4class{GenomicRanges}} or 
    \code{\linkS4class{GRangesList}} object used to map between \code{x} and 
    the result. The ranges can be any feature in the transcriptome extracted
    from a \code{TxDb} (e.g., introns, exons, cds regions). See 
    ?\code{transcripts} and ?\code{transcriptsBy} for a list of extractor 
    functions.
  }
  \item{reverse}{
    TRUE or FALSE. When FALSE (default) the range in \code{x} is mapped to the
    transcript space. Seqlevels of the output object are the names on the
    \code{transcripts} object; most often these will be the transcript names. 
    
    When \code{reverse=TRUE}, \code{x} is mapped to genome space. Seqlevels
    of the output are the chromosome names.
  }
  \item{ignore.strand}{
    When TRUE, strand is ignored in overlap operations. When
    \code{transcripts} is a \code{GRangesList} and \code{ignore.strand =
    FALSE} all inner list elements of a common list element must have the same
    strand. 
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
  \itemize{
    \item{mapToTranscripts:}{
      \itemize{
        \item{}{
          When \code{reverse=FALSE} (default), the genomic range in \code{x} is 
          mapped to the local position in the \code{transcripts} ranges. A 
          successful mapping occurs when \code{x} is completely within the 
          \code{transcripts} range, equivalent to: 
          \preformatted{
          findOverlaps(..., type="within")
          }
          Transcriptome-based coordinates start counting at 1 at the beginning
          of the \code{transcripts} range and return positions where \code{x}
          was aligned. The seqlevels of the return object are taken from the
          \code{transcripts} object and should be transcript names. In this
          direction, mapping is attempted between all elements of \code{x} and
          all elements of \code{transcripts}.
        }
        \item{}{
          When \code{reverse=TRUE}, the transcripts range in \code{x} is 
          mapped to genomic coordinates using the ranges in \code{transcripts}.
          A successful mapping occurs when the following is TRUE:
          \preformatted{
          width(transcripts) >= start(x) + width(x)
          }
          \code{x} is aligned to \code{transcripts} by moving in
          \code{start(x)} positions in from the beginning of the
          \code{transcripts} range.  The seqlevels of the return object are
          genomic chromosomes.

          When mapping to the genome, name matching is used to determine the
          mapping pairs (vs attempting to match all possible pairs). Ranges in
          \code{x} are only mapped to ranges in \code{transcripts} with the
          same name. Name matching is motivated by use cases such as
          differentially expressed regions where the expressed regions in
          \code{x} would only be related to a subset of regions in
          \code{transcripts}, which may contains gene or transcript ranges.
        }
      }
    }
    \item{pmapToTranscripts:}{
      Element-wise (aka `parallel`) version of \code{mapToTranscripts}. The
      i-th range in \code{x} is mapped to the i-th range in
      \code{transcripts}; \code{x} and \code{transcripts} must have the same
      length. When \code{reverse=FALSE} mapping is to the transcripts;
      when TRUE mapping is to the genome.
 
      Ranges in \code{x} that do not map (out of bounds or strand mismatch) 
      are returned as zero-width ranges starting at 1. These ranges are given
      the special seqname of "unmapped". Note the non-parallel methods do not 
      return unmapped ranges so the "unmapped" seqname is unique to 
      \code{pmapToTranscripts}.
    }
  }
}

\value{
  An object the same class as \code{x}.

  Parallel methods return an object the same shape as \code{x}. Ranges that
  are non-hits (out of bounds or strand-mismatch) are returned as zero-width
  ranges starting at 1. The seqname of the non-hits is "unmapped".

  Non-parallel methods return an object that varies in length similar to a
  Hits object. The result only contains mapped records, strand mismatch
  and non-hits are not returned. \code{xHits} and \code{transcriptsHits}
  metadata columns record the element of \code{x} and \code{transcripts} used
  in the mapping.

  When \code{reverse=FALSE} seqlevels of the output object are the names on 
  the \code{transcripts} object; most often these will be transcript names. 
  When \code{reverse=TRUE} seqlevels of the output are chromosome names.
}

\seealso{
  \itemize{
    \item \code{?\link[GenomicAlignments]{mapToAlignments}} in the 
    \pkg{GenomicAlignments} package for methods mapping between
    reads and genome space using a CIGAR alignment.
  }
}

\author{V. Obenchain, M. Lawrence and H. Pages}

\examples{

## ---------------------------------------------------------------------
## A. Basic Use
## ---------------------------------------------------------------------

library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
cds <- cdsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, "tx", use.names=TRUE)[1:3]
x <- GRanges("chr2L", 
             IRanges(c(7500, 8400, 9000), width=200, names=LETTERS[1:3])) 

## Map to transcript coordinates:
mapToTranscripts(x, cds)

## Note the seqnames of the output are the transcript names, not
## chromosomes. The 'xHits' and 'transcriptsHits' metadata 
## columns indicate which elements were involved in the mapping.

## The element-wise version returns both mapped and unmapped ranges.
x <- GRanges("chr2L", IRanges(c(8100, 8600), width=10, names=LETTERS[1:2]))
pmapToTranscripts(x, cds[1:2])

## Mapping to genome space requires both 'x' and 'transcripts' to have
## names. A map is only attempted for ranges with matching names.
x <- GRanges("chr1", 
             IRanges(c(1, 20, 10, 3), width=5, names=c("A", "B", "A", "C")))
gr <- GRanges("chr1", 
              IRanges(c(25, 30, 3), width=10, names=c("C", "C", "A"))) 
mapToTranscripts(x, gr, reverse = TRUE)

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
pmapToTranscripts(first_loc, exons[[1]], reverse = TRUE)


## -----------------------------------------------------------------------
## C. Map 3'UTR variants to genome coordinates
## -----------------------------------------------------------------------

## A study by Skeeles et. al (PLoS ONE 8(3): e58609. doi:
## 10.1371/journal.pone.0058609) investigated the impact of 3'UTR variants 
## on the expression of cancer susceptibility genes.

## 8 candidate miRNA genes on chromosome 12 were used to test for 
## differential luciferase expression in mice. In Table 2 of the manuscript
## variant locations are given as nucleotide position within the gene.
geneNames <- c("Bcap29", "Dgkb", "Etv1", "Hbp1", "Hbp1", "Ifrd1", "Ifrd1", 
               "Pik3cg", "Pik3cg", "Tspan13", "Twistnb")
starts <- c(1409, 3170, 3132, 2437, 2626, 3239, 3261, 4947, 4979, 958, 1489) 
snps <- GRanges("chr12", IRanges(starts, width=1, names=geneNames))

## To map these transcript-space coordinates to the genome we need gene ranges
## in genome space.
library(org.Mm.eg.db)
geneid <- select(org.Mm.eg.db, unique(geneNames), "ENTREZID", "SYMBOL")
geneid

## Extract the gene regions:
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)[geneid$ENTREZID]

## mapToTranscripts(..., reverse=TRUE) determines pairs to map by comparing 
## names in 'x' and 'transcripts'. Ranges in 'snps' will be mapped to all 
## ranges in 'genes' with the same name. Currently the names of 'genes' are 
## internal gene ids. Rename 'genes' with the appropriate gene symbol.
names(genes) <- geneid$SYMBOL 

## The xHits and transcriptsHits metadta columns indicate which ranges in
## 'snps' and 'genes' were involved in the mapping.
mapToTranscripts(snps, genes, reverse = TRUE)


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
mapToTranscripts(snps, unlist(cDNA), reverse = FALSE)

## Map the first SNP to transcript uc002pct.1 and the second to uc002pcu.1.
pmapToTranscripts(snps[1:2], unlist(cDNA)[1:2], reverse = FALSE)

## The cdsBy() extractor returns coding regions by gene or by transcript.
## Extract the coding regions for transcript uc002pct.1.
cds <- cdsBy(txdb, "tx", use.names=TRUE)["uc002pct.1"]
cds

## The 'cds' object is a GRangesList of length 1 containing all cds ranges
## for the single transcript uc002pct.1. To map to the concatenated group
## of ranges we leave 'cds' as a GRangesList.

## Map all 5 SNPs to the total cds region: 
mapToTranscripts(snps, cds, reverse = FALSE)

## Only the second SNP could be mapped. Unlisting the 'cds' object maps the
## SNPs to the individual cds ranges (vs the concatenated range).
mapToTranscripts(snps[2], unlist(cds), reverse = FALSE)

## The location is the same because the SNP hit the first cds range. If the
## transcript had been on the negative strand the difference in mapping

## would be more obvious.
strand(snps) <- strand(cds) <- "-"
mapToTranscripts(snps[2], cds, reverse = FALSE, ignore.strand=FALSE)
mapToTranscripts(snps[2], unlist(cds), reverse = FALSE, ignore.strand=FALSE)
}

\keyword{methods}
\keyword{utilities}
