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
                        chrominfo,
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
  \item{chrominfo}{data frame containing information about the
  chromosomes.  This data.frame has 3 columns: 'chrom', 'length' and
  'is_circular'.  The 'chrom' column contains the character names of the
  chromosome elements that are present, the 'length' is the integer
  length for each one and the 'is_circular' contains a logical value to
  indicate whether that element is a circularized (such as a
  mitochondria or a chloroplast).  If this argument is left blank then
  there will not be any length information recorded for the different
  chromosomes since it is not possible to infer if from the transcript
  ranges alone.}
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
  information.  For GTF files this is usually given as "exon_number",
  however you still must specify this argument if you don't want the
  code to try and infer the exon rank information.  For gff3 files, we
  have not seen any examples of this information encoded anywhere, but
  if you have a file with an attribute, you can still specify this to
  avoid the inference.
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
         circ_seqs=DEFAULT_CIRC_SEQS,
         miRBaseBuild=NULL)
if(interactive()) {
saveDb(txdb,file="TESTGFF.sqlite")
}


## TESTING GTF this time specifying the chrominfo
gtfFile <- system.file("extdata","Aedes_aegypti.partial.gtf",
                       package="GenomicFeatures")
chrominfo <- data.frame(chrom = c('supercont1.1','supercont1.2'),
                        length=c(5220442, 5300000),
                        is_circular=c(FALSE, FALSE))
txdb2 <- makeTranscriptDbFromGFF(file=gtfFile,
         format="gtf",
         gffExonRankAttributeName="exon_number",
         chrominfo= chrominfo,
         dataSource=paste("ftp://ftp.ensemblgenomes.org/pub/metazoa/",
                          "release-13/gtf/aedes_aegypti/",sep=""),
         species="Aedes aegypti",
         circ_seqs=DEFAULT_CIRC_SEQS,
         miRBaseBuild=NULL)
if(interactive()) {
saveDb(txdb2,file="TESTGTF.sqlite")
}

}
