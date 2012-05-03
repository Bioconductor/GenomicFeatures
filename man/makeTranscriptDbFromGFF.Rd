\name{makeTranscriptDbFromGFF}

\alias{makeTranscriptDbFromGFF}

\title{
  Make a TranscriptDb object from annotations available as GFF3 or GTF
  formatted files
}
\description{
  The \code{makeTranscriptDbFromGFF} function allows the user
  to make a \link{TranscriptDb} object from transcript annotations
  available as GFF3 or GTF files.
}
\usage{
makeTranscriptDbFromGFF(file,
                        format=c("gff3","gtf"),
                        gffExonRankAttributeName,
                        dataSource,
                        species,
                        circ_seqs=DEFAULT_CIRC_SEQS,
                        miRBaseBuild = NULL)
}
\arguments{
  \item{file}{path/file to be processed}
  \item{format}{"gff3" or "gtf" depending on which file format you have
  to process}
  \item{gffExonRankAttributeName}{what is the name (if any) of the
  attribute that defines the exon rank information.}
  \item{dataSource}{Where did this data file originate? Please be as
  specific as possible.}
  \item{species}{What is the Genus and species of this organism. Please
  use proper scientific nomenclature for example: "Homo sapiens" or
  "Canis familiaris" and not "human" or "fuzzy buddy".  If properly
  written, this information may be used by the software to help you out
  later.}
  \item{circ_seqs}{a character vector to list out which chromosomes
    should be marked as circular.}
  \item{miRBaseBuild}{specify the string for the appropriate build
    Information from mirbase.db to use for microRNAs.  This can be
    learned by calling \code{supportedMiRBaseBuildValues}.  By default,
    this value will be NULL, which will inactivate the \code{microRNAs}
    accessor.}
}
\details{
  \code{makeTranscriptDbFromGFF} is a convenience function that feeds
  data from the parsed file to the lower level \code{\link{makeTranscriptDb}}
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
  information.  For GTF files this is usually given as "exon_number" and
  so if you don't specify this and are processing a gtf file, then we
  will assume that the exon rank information is stored in an attribute
  called "exon_number".  For gff3 files, we have not seen any examples
  of this information encoded anywhere and so we will normally try to infer
  this information from the position along the chromosome.  However, for
  gff3 files, if you have a file that actually does have an attribute
  defined for the exon rank, you can use the gffExonRankAttributeName to
  specify it.
}

\value{A \link{TranscriptDb} object.}

\author{
  M. Carlson and H. Pages
}

\seealso{
  \code{\link{DEFAULT_CIRC_SEQS}},
  \code{\link{makeTranscriptDb}}
  \code{\link{supportedMiRBaseBuildValues}}
}

\examples{

## TESTING GFF3
gffFile <- system.file("extdata","a.gff3",package="GenomicFeatures")
txdb <- makeTranscriptDbFromGFF(file=gffFile,
         format="gff3",
         dataSource="partial gtf file for Tomatoes for testing",
         species="Solanum lycopersicum",
         gffExonRankAttributeName = "nb_exon",
         circ_seqs=DEFAULT_CIRC_SEQS,
         miRBaseBuild=NULL)
if(interactive()) {
saveDb(txdb,file="TESTGFF.sqlite")
}


## TESTING GTF
gtfFile <- system.file("extdata","Aedes_aegypti.partial.gtf",
                       package="GenomicFeatures")
txdb <- makeTranscriptDbFromGFF(file=gtfFile,
         format="gtf",
         dataSource=paste("ftp://ftp.ensemblgenomes.org/pub/metazoa/",
                          "release-13/gtf/aedes_aegypti/",sep=""),
         species="Aedes aegypti",
         circ_seqs=DEFAULT_CIRC_SEQS,
         miRBaseBuild=NULL)
if(interactive()) {
saveDb(txdb,file="TESTGTF.sqlite")
}

}
