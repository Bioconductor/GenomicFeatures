\name{makeTxDbPackage}

\alias{makeTxDbPackage}
\alias{makeTxDbPackageFromUCSC}
\alias{makeTxDbPackageFromBiomart}


\title{
  Making a TranscriptDb packages from annotations available at the
  UCSC Genome Browser, biomaRt or from another source.
}
\description{
  The \code{makeTxDbPackageFromUCSC} function allows the user
  to make a \link{TranscriptDb} object from transcript annotations
  available at the UCSC Genome Browser.
  The \code{makeTxDbPackageFromBiomart} function allows the user
  to do the same thing as \code{makeTxDbPackageFromUCSC} except that the
  annotations originate from biomaRt.
  Finally, the \code{makeTxDbPackage} function allows the user to make a   
  \link{TranscriptDb} object from transcript annotations that are in 
  a custom transcript Database, such as could be produced using 
  \code{makeTranscriptDb}.
}
\usage{

  makeTxDbPackageFromUCSC(
    version=,
    maintainer,
    author,
    destDir=".",
    license="Artistic-2.0",
    genome="hg19",
    tablename="knownGene",
    transcript_ids=NULL,
    circ_seqs=DEFAULT_CIRC_SEQS,
    url="http://genome.ucsc.edu/cgi-bin/",
    goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath")
  
  makeTxDbPackageFromBiomart(
    version,
    maintainer,
    author,
    destDir=".",
    license="Artistic-2.0",
    biomart="ensembl",
    dataset="hsapiens_gene_ensembl",
    transcript_ids=NULL,
    circ_seqs=DEFAULT_CIRC_SEQS)

  makeTxDbPackage(txdb,
                  version,
  		            maintainer,
                  author,
  	              destDir=".",
                  license="Artistic-2.0")
}
\arguments{
  \item{version}{What is the version number for this package?}
  \item{maintainer}{Who is the package maintainer? (must include email to be 
    valid)}
  \item{author}{Who is the creator of this package?}
  \item{destDir}{A path where the package source should be assembled.}
  \item{license}{What is the license (and it's version)}
  \item{biomart}{which BioMart database to use.
    Get the list of all available BioMart databases with the
    \code{\link[biomaRt]{listMarts}} function from the biomaRt
    package. See the details section below for a list of BioMart
    databases with compatible transcript annotations.}
  \item{dataset}{which dataset from BioMart. For example:
    \code{"hsapiens_gene_ensembl"}, \code{"mmusculus_gene_ensembl"},
    \code{"dmelanogaster_gene_ensembl"}, \code{"celegans_gene_ensembl"},
    \code{"scerevisiae_gene_ensembl"}, etc in the ensembl database.
    See the examples section below for how to discover which datasets
    are available in a given BioMart database.}
  \item{genome}{genome abbreviation used by UCSC and obtained by
    \code{\link[rtracklayer]{ucscGenomes}()[ , "db"]}.
    For example: \code{"hg18"}.}
  \item{tablename}{name of the UCSC table containing the transcript
    annotations to retrieve. Use the \code{supportedUCSCtables} utility
    function to get the list of supported tables. Note that not all
    tables are available for all genomes.}
  \item{transcript_ids}{optionally, only retrieve transcript
    annotation data for the specified set of transcript ids.
    If this is used, then the meta information displayed for the
    resulting \link{TranscriptDb} object will say 'Full dataset: no'.
    Otherwise it will say 'Full dataset: yes'.}
  \item{circ_seqs}{a character vector to list out which chromosomes
    should be marked as circular.}
  \item{url,goldenPath_url}{use to specify the location of an
    alternate UCSC Genome Browser.}
  \item{txdb}{A \link{TranscriptDb} object that represents a handle to a 
    transcript database. This object type is what is returned by 
    \code{makeTranscriptDbFromUCSC}, \code{makeTranscriptDbFromUCSC} or 
    \code{makeTranscriptDb}}
}
\details{
  \code{makeTxDbPackageFromUCSC} is a convenience function that calls
    both the \code{\link{makeTranscriptDbFromUCSC}} and the             
    \code{\link{makeTxDbPackage}} functions.  The 
    \code{makeTxDbPackageFromBiomart} follows a similar pattern and calls the       
    \code{\link{makeTranscriptDbFromBiomart}} and \code{\link{makeTxDbPackage}} 
    functions.  
}

\value{A \link{TranscriptDb} object.}

\author{
  M. Carlson
}

\seealso{
  \code{\link[rtracklayer]{ucscGenomes}},
  \code{\link{DEFAULT_CIRC_SEQS}},
  \code{\link{makeTranscriptDbFromUCSC}},
  \code{\link{makeTranscriptDbFromBiomart}},
  \code{\link{makeTranscriptDb}}
}

\examples{

## Display the list of tables supported by makeTxDbPackageFromUCSC():
supportedUCSCtables()

\dontrun{
## Makes a transcript package for Yeast from the ensGene table at UCSC:
makeTxDbPackageFromUCSC(version="0.01", 
                        maintainer="Some One <so@someplace.org>", 
                        author="Some One <so@someplace.com>",
                        genome="sacCer2", 
                        tablename="ensGene")

## Makes a transcript package from Human by using biomaRt and limited to a 
## small subset of the transcripts.
transcript_ids <- c(
    "ENST00000400839",
    "ENST00000400840",
    "ENST00000478783",
    "ENST00000435657",
    "ENST00000268655",
    "ENST00000313243",
    "ENST00000341724")
    
makeTxDbPackageFromBiomart(version="0.01", 
                           maintainer="Some One <so@someplace.org>", 
                           author="Some One <so@someplace.com>",
                           transcript_ids=transcript_ids)

}

}
