\name{makeTxDbPackage}

\alias{makeTxDbPackage}
\alias{makeTxDbPackageFromUCSC}
\alias{makeFDbPackageFromUCSC}
\alias{makeTxDbPackageFromBiomart}
\alias{supportedMiRBaseBuildValues}

\title{
  Making a TxDb package from annotations available at the
  UCSC Genome Browser, biomaRt or from another source.
}
\description{
  A \link{TxDb} package is an annotation package containing a
  \link{TxDb} object.

  The \code{makeTxDbPackageFromUCSC} function allows the user
  to make a \link{TxDb} package from transcript annotations
  available at the UCSC Genome Browser.

  The \code{makeTxDbPackageFromBiomart} function allows the user
  to do the same thing as \code{makeTxDbPackageFromUCSC} except that the
  annotations originate from biomaRt.

  Finally, the \code{makeTxDbPackage} function allows the user to make a   
  \link{TxDb} package directly from a \link{TxDb} object.
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
    goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
    taxonomyId=NA,
    miRBaseBuild=NA)

makeFDbPackageFromUCSC(
    version,
    maintainer,
    author,
    destDir=".",
    license="Artistic-2.0",
    genome="hg19",
    track="tRNAs",
    tablename="tRNAs",
    columns = UCSCFeatureDbTableSchema(genome, track, tablename),
    url="http://genome.ucsc.edu/cgi-bin/",
    goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
    chromCol=NULL,
    chromStartCol=NULL,
    chromEndCol=NULL,
    taxonomyId=NA)

  makeTxDbPackageFromBiomart(
    version,
    maintainer,
    author,
    destDir=".",
    license="Artistic-2.0",
    biomart="ENSEMBL_MART_ENSEMBL",
    dataset="hsapiens_gene_ensembl",
    transcript_ids=NULL,
    circ_seqs=DEFAULT_CIRC_SEQS,
    filter=NULL,
    id_prefix="ensembl_",
    host="www.ensembl.org",
    port=80,
    taxonomyId=NA,
    miRBaseBuild=NA)

  makeTxDbPackage(txdb,
                  version,
  		  maintainer,
                  author,
  	          destDir=".",
                  license="Artistic-2.0",
                  pkgname=NULL)

  supportedMiRBaseBuildValues()
}
\arguments{
  \item{version}{What is the version number for this package?}
  \item{maintainer}{Who is the package maintainer? (must include email
    to be valid). Should be a \code{\link{person}} object, or something
    coercible to one, like a string. May be omitted if the \code{author}
    argument is a \code{person} containing someone with the maintainer role.}
  \item{author}{Who is the creator of this package? Should be
    a \code{\link{person}} object, or something coercible to one, like a
    character vector of names. The \code{maintainer} argument will be
    merged into this list.}
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
  \item{track}{name of the UCSC track.  Use
    \code{supportedUCSCFeatureDbTracks} to get the list of available
        tracks for a particular genome}
  \item{tablename}{name of the UCSC table containing the transcript
    annotations to retrieve. Use the \code{supportedUCSCtables} utility
    function to get the list of supported tables. Note that not all
    tables are available for all genomes.}
  \item{transcript_ids}{optionally, only retrieve transcript
    annotation data for the specified set of transcript ids.
    If this is used, then the meta information displayed for the
    resulting \link{TxDb} object will say 'Full dataset: no'.
    Otherwise it will say 'Full dataset: yes'.}
  \item{circ_seqs}{a character vector to list out which chromosomes
    should be marked as circular.}
  \item{filter}{Additional filters to use in the BioMart query. Must be
    a named list. An example is \code{filter=as.list(c(source="entrez"))}}
  \item{host}{The host URL of the BioMart. Defaults to www.ensembl.org.}
  \item{port}{The port to use in the HTTP communication with the host.}
  \item{id_prefix}{Specifies the prefix used in BioMart attributes. For
    example, some BioMarts may have an attribute specified as
    \code{"ensembl_transcript_id"} whereas others have the same attribute
    specified as \code{"transcript_id"}. Defaults to \code{"ensembl_"}.}
  \item{columns}{a named character vector to list out the names and
        types of the other columns that the downloaded track should
        have.  Use \code{UCSCFeatureDbTableSchema} to retrieve this
        information for a particular table.}
  \item{url,goldenPath_url}{use to specify the location of an
    alternate UCSC Genome Browser.}
  \item{chromCol}{If the schema comes back and the 'chrom' column has been
        labeled something other than 'chrom', use this argument to indicate
        what that column has been labeled as so we can properly designate it.
        This could happen (for example) with the knownGene track tables, which
        has no 'chromStart' or 'chromEnd' columns, but which DOES have columns
        that could reasonably substitute for these columns under particular
        circumstances.  Therefore we allow these three columns to have arguments
        so that their definition can be re-specified}
  \item{chromStartCol}{Same thing as chromCol, but for renames of 'chromStart'}
  \item{chromEndCol}{Same thing as chromCol, but for renames of 'chromEnd'}
  \item{txdb}{A \link{TxDb} object that represents a handle to a 
    transcript database. This object type is what is returned by 
    \code{makeTxDbFromUCSC}, \code{makeTxDbFromUCSC} or 
    \code{makeTxDb}}
  \item{taxonomyId}{By default this value is NA and the organism
    provided (or inferred) will be used to look up the correct value for
    this.  But you can use this argument to override that and supply
    your own valid taxId here}
  \item{miRBaseBuild}{specify the string for the appropriate build
    Information from mirbase.db to use for microRNAs.  This can be
    learned by calling \code{supportedMiRBaseBuildValues}.  By default,
    this value will be set to \code{NA}, which will inactivate the
    \code{microRNAs} accessor.}
  \item{pkgname}{By default this value is NULL and does not need to be
    filled in (a package name will be generated for you).  But if you
    override this value, then the package and it's object will be
    instead named after this value.  Be aware that the standard rules
    for package names will apply, (so don't include spaces,
    underscores or dashes)}  
}
\details{
  \code{makeTxDbPackageFromUCSC} is a convenience function that calls
    both the \code{\link{makeTxDbFromUCSC}} and the
    \code{\link{makeTxDbPackage}} functions.  The
    \code{makeTxDbPackageFromBiomart} follows a similar pattern and
    calls the \code{\link{makeTxDbFromBiomart}} and
    \code{\link{makeTxDbPackage}} functions.
  \code{supportedMiRBaseBuildValues} is a convenience function that will
  list all the possible values for the miRBaseBuild argument.
}

\value{A \link{TxDb} object.}

\author{
  M. Carlson
}

\seealso{
  \code{\link{makeTxDbFromUCSC}},
  \code{\link{makeTxDbFromBiomart}},
  \code{\link{makeTxDb}},
  \code{\link[rtracklayer]{ucscGenomes}},
  \code{\link{DEFAULT_CIRC_SEQS}}
}

\examples{
## First consider relevant helper/discovery functions:
## Display the list of tables supported by makeTxDbPackageFromUCSC():
supportedUCSCtables()

## Can also list all the possible values for the miRBaseBuild argument:
supportedMiRBaseBuildValues()

## Next are examples of actually building a package:
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
