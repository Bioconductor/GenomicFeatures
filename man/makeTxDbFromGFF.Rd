\name{makeTxDbFromGFF}

\alias{makeTxDbFromGFF}

% deprecated
\alias{makeTranscriptDbFromGFF}

\title{
  Make a TxDb object from annotations available as a
  GFF3 or GTF file
}
\description{
  The \code{makeTxDbFromGFF} function allows the user
  to make a \link{TxDb} object from transcript annotations
  available as a GFF3 or GTF file.
}
\usage{
makeTxDbFromGFF(file,
                format=c("gff3","gtf"),
                exonRankAttributeName=NA,
                gffGeneIdAttributeName=NA,
                chrominfo=NA,
                dataSource=NA,
                organism=NA,
                circ_seqs=DEFAULT_CIRC_SEQS,
                miRBaseBuild=NA,
                useGenesAsTranscripts=FALSE,
                gffTxName="mRNA",
                species=NA)
}
\arguments{
  \item{file}{path/file to be processed}
  \item{format}{"gff3" or "gtf" depending on which file format you have
  to process}
  \item{exonRankAttributeName}{ignored and deprecated}
  \item{gffGeneIdAttributeName}{ignored and deprecated}
  \item{chrominfo}{data frame containing information about the
  chromosomes.  Will be passed to the internal call to
  \code{\link{makeTxDb}}. See \code{?\link{makeTxDb}} for the details.}
  \item{dataSource}{Where did this data file originate?  Please be as
  specific as possible.}
  \item{organism}{What is the Genus and species of this organism. Please
  use proper scientific nomenclature for example: "Homo sapiens" or
  "Canis familiaris" and not "human" or "my fuzzy buddy".  If properly
  written, this information may be used by the software to help you out
  later.}
  \item{circ_seqs}{a character vector to list out which chromosomes
    should be marked as circular.}
  \item{miRBaseBuild}{specify the string for the appropriate build
    Information from mirbase.db to use for microRNAs.  This can be
    learned by calling \code{supportedMiRBaseBuildValues}.  By default,
    this value will be set to \code{NA}, which will inactivate the
    \code{microRNAs} accessor.}
  \item{useGenesAsTranscripts}{ignored and deprecated}
  \item{gffTxName}{ignored and deprecated}
  \item{species}{deprecated in favor of \code{organism}}
}
\details{
  \code{makeTxDbFromGFF} is a convenience function that feeds
  data from the parsed file to the \code{\link{makeTxDbFromGRanges}}
  function.
}

\value{A \link{TxDb} object.}

\author{
  M. Carlson and H. Pages
}

\seealso{
  \itemize{
    \item \code{\link{makeTxDbFromGRanges}}, which \code{makeTxDbFromGFF}
          is based on, for making a \link{TxDb} object from a
          \link[GenomicRanges]{GRanges} object.

    \item The \code{\link[rtracklayer]{import}} function in the
          \pkg{rtracklayer} package (also used by \code{makeTxDbFromGFF}
          internally).

    \item \code{\link{makeTxDbFromUCSC}} and \code{\link{makeTxDbFromBiomart}}
          for convenient ways to make a \link{TxDb} object from UCSC or BioMart
          online resources.

    \item \code{\link{DEFAULT_CIRC_SEQS}}.

    \item The \code{\link{supportedMiRBaseBuildValues}} function for
          listing all the possible values for the \code{miRBaseBuild}
          argument.

    \item The \link{TxDb} class.

    \item \code{\link{makeTxDb}} for the low-level function used by the
          \code{makeTxDbFrom*} functions to make the \link{TxDb} object
          returned to the user.
  }
}

\examples{
## TESTING GFF3
gffFile <- system.file("extdata","GFF3_files","a.gff3",package="GenomicFeatures")
txdb <- makeTxDbFromGFF(file=gffFile,
            format="gff3",
            dataSource="partial gtf file for Tomatoes for testing",
            organism="Solanum lycopersicum")

## TESTING GTF, this time specifying the chrominfo
gtfFile <- system.file("extdata","GTF_files","Aedes_aegypti.partial.gtf",
                       package="GenomicFeatures")
chrominfo <- data.frame(chrom = c('supercont1.1','supercont1.2'),
                        length=c(5220442, 5300000),
                        is_circular=c(FALSE, FALSE))
txdb2 <- makeTxDbFromGFF(file=gtfFile,
             format="gtf",
             chrominfo=chrominfo,
             dataSource=paste("ftp://ftp.ensemblgenomes.org/pub/metazoa/",
                              "release-13/gtf/aedes_aegypti/",sep=""),
             organism="Aedes aegypti")
}
