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
                species=NA,
                circ_seqs=DEFAULT_CIRC_SEQS,
                miRBaseBuild=NA,
                useGenesAsTranscripts=FALSE,
                gffTxName="mRNA")
}
\arguments{
  \item{file}{path/file to be processed}
  \item{format}{"gff3" or "gtf" depending on which file format you have
  to process}
  \item{exonRankAttributeName}{character(1) name of the attribute that
    defines the exon rank information, or \code{NA} to indicate that
    exon ranks are inferred from order of occurrence in the GFF.}
  \item{gffGeneIdAttributeName}{an optional argument that can be used
  for gff style files ONLY.  If the gff file lacks rows to specify gene
  IDs but the mRNA rows of the gff file specify the gene IDs via a named
  attribute,then passing the name of the attribute for this argument can
  allow the file to still extract gene IDs that map to these
  transcripts.  If left blank, then the parser will try and extract rows
  that are named 'gene' for gene to transcript mappings when parsing a
  gff3 file.  For gtf files this argument is ignored entirely.}
  \item{chrominfo}{data frame containing information about the
  chromosomes.  Will be passed to the internal call to
  \code{\link{makeTxDb}}. See \code{?\link{makeTxDb}} for the details.}
  \item{dataSource}{Where did this data file originate?  Please be as
  specific as possible.}
  \item{species}{What is the Genus and species of this organism. Please
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
  \item{useGenesAsTranscripts}{This flag is normally off, but if enabled
    it will try to salvage a file that has no RNA features by assuming
    that you can use the ranges available for the Gene features in their
    place.  Obviously, this is something you won't want to do unless you
    are dealing with something very simple like a prokaryote.}
  \item{gffTxName}{If you are parsing a gff file and the transcript
    feature type is named after something other than 'mRNA' (the default),
    then you can specify an alternate string value to parse with here.}
}
\details{
  \code{makeTxDbFromGFF} is a convenience function that feeds
  data from the parsed file to the lower level \code{\link{makeTxDb}}
  function.

  There are some real deficiencies in the gtf and the gff3 file formats
  to bear in mind when making use of them.  For gtf files the length of
  the transcripts is not normally encoded and so it has to be inferred
  from the exon ranges presented.  That's not a horrible problem, but it
  bears mentioning for the sake of full disclosure.  And for gff3 files
  the situation is typically even worse since they usually don't encode
  any information about the exon rank within a transcript.  This is a
  serious oversight and so if you have an alternative to using this kind
  of data, you should really do so.
  
  Some files will have an attribute defined to indicate the exon rank
  information.  For GTF files this is usually given as "exon_number",
  however you still must specify this argument if you don't want the
  code to try and infer the exon rank information.  For gff3 files, we
  have not seen any examples of this information encoded anywhere, but
  if you have a file with an attribute, you can still specify this to
  avoid the inference.
}

\value{A \link{TxDb} object.}

\author{
  M. Carlson
}

\seealso{
  \code{\link{DEFAULT_CIRC_SEQS}},
  \code{\link{makeTxDbFromUCSC}},
  \code{\link{makeTxDbFromBiomart}},
  \code{\link{makeTxDb}},
  \code{\link{supportedMiRBaseBuildValues}}
}

\examples{
## TESTING GFF3
gffFile <- system.file("extdata","a.gff3",package="GenomicFeatures")
txdb <- makeTxDbFromGFF(file=gffFile,
            format="gff3",
            exonRankAttributeName=NA,
            dataSource="partial gtf file for Tomatoes for testing",
            species="Solanum lycopersicum")
if(interactive()) {
saveDb(txdb,file="TESTGFF.sqlite")
}

## TESTING GTF, this time specifying the chrominfo
gtfFile <- system.file("extdata","Aedes_aegypti.partial.gtf",
                       package="GenomicFeatures")
chrominfo <- data.frame(chrom = c('supercont1.1','supercont1.2'),
                        length=c(5220442, 5300000),
                        is_circular=c(FALSE, FALSE))
txdb2 <- makeTxDbFromGFF(file=gtfFile,
             format="gtf",
             exonRankAttributeName="exon_number",
             chrominfo=chrominfo,
             dataSource=paste("ftp://ftp.ensemblgenomes.org/pub/metazoa/",
                              "release-13/gtf/aedes_aegypti/",sep=""),
             species="Aedes aegypti")
if(interactive()) {
    saveDb(txdb2,file="TESTGTF.sqlite")
}
}
