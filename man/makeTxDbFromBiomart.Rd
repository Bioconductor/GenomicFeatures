\name{makeTxDbFromBiomart}

\alias{makeTxDbFromBiomart}
\alias{getChromInfoFromBiomart}

\title{
  Make a TxDb object from annotations available on a
  BioMart database
}
\description{
  The \code{makeTxDbFromBiomart} function allows the user
  to make a \link{TxDb} object from transcript annotations
  available on a BioMart database.
}
\usage{
makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                    dataset="hsapiens_gene_ensembl",
                    transcript_ids=NULL,
                    circ_seqs=DEFAULT_CIRC_SEQS,
                    filters="",
                    id_prefix="ensembl_",
                    host="www.ensembl.org",
                    port=80,
                    taxonomyId=NA,
                    miRBaseBuild=NA)

getChromInfoFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                        dataset="hsapiens_gene_ensembl",
                        id_prefix="ensembl_",
                        host="www.ensembl.org",
                        port=80)
}
\arguments{
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
  \item{transcript_ids}{optionally, only retrieve transcript
    annotation data for the specified set of transcript ids.
    If this is used, then the meta information displayed for the
    resulting \link{TxDb} object will say 'Full dataset: no'.
    Otherwise it will say 'Full dataset: yes'.}
  \item{circ_seqs}{a character vector to list out which chromosomes
    should be marked as circular.}
  \item{filters}{Additional filters to use in the BioMart query. Must be
    a named list. An example is \code{filters=list(source="entrez")}}
  \item{id_prefix}{Specifies the prefix used in BioMart attributes. For
    example, some BioMarts may have an attribute specified as
    \code{"ensembl_transcript_id"} whereas others have the same attribute
    specified as \code{"transcript_id"}. Defaults to \code{"ensembl_"}.}
  \item{host}{The host URL of the BioMart. Defaults to www.ensembl.org.}
  \item{port}{The port to use in the HTTP communication with the host.}
  \item{taxonomyId}{By default this value is NA and the dataset
    selected will be used to look up the correct value for this.  But
    you can use this argument to override that and supply your own
    taxId here (which will be independently checked to make sure its a
    real taxonomy id).  Normally you should never need to use this.}  
  \item{miRBaseBuild}{specify the string for the appropriate build
    Information from mirbase.db to use for microRNAs.  This can be
    learned by calling \code{supportedMiRBaseBuildValues}.  By default,
    this value will be set to \code{NA}, which will inactivate the
    \code{microRNAs} accessor.}
}

\details{
  \code{makeTxDbFromBiomart} is a convenience function that feeds
  data from a BioMart database to the lower level
  \code{\link{makeTxDb}} function.
  See \code{?\link{makeTxDbFromUCSC}} for a similar function
  that feeds data from the UCSC source.

  Here is a list of datasets known to be compatible with
  \code{makeTxDbFromBiomart} (updated on November 12, 2015):
  \itemize{
    \item All the datasets in the main Ensembl database (list them with
          \code{biomaRt::listDatasets(biomaRt::useMart(
          biomart="ENSEMBL_MART_ENSEMBL", host="www.ensembl.org"))}).

    \item All the datasets in the Ensembl Fungi database (list them with
          \code{biomaRt::listDatasets(biomaRt::useMart(
          biomart="fungal_mart", host="fungi.ensembl.org"))}).

    \item All the datasets in the Ensembl Metazoa database (list them with
          \code{biomaRt::listDatasets(biomaRt::useMart(
          biomart="metazoa_mart", host="metazoa.ensembl.org"))}).

    \item All the datasets in the Ensembl Plants database (list them with
          \code{biomaRt::listDatasets(biomaRt::useMart(
          biomart="plants_mart", host="plants.ensembl.org"))}).

    \item All the datasets in the Ensembl Protists database (list them with
          \code{biomaRt::listDatasets(biomaRt::useMart(
          biomart="protist_mart", host="protists.ensembl.org"))}).

    \item All the datasets in the Gramene Mart (list them with
          \code{biomaRt::listDatasets(biomaRt::useMart(
          biomart="ENSEMBL_MART_PLANT", host="ensembl.gramene.org"))}).
  }
  Not all these datasets have CDS information.
}

\value{
  A \link{TxDb} object for \code{makeTxDbFromBiomart}.

  A data frame with 1 row per chromosome (or scaffold) and with columns
  \code{chrom} and \code{length} for \code{getChromInfoFromBiomart}.
}

\author{
  M. Carlson and H. Pagès
}

\seealso{
  \itemize{
    \item \code{\link{makeTxDbFromUCSC}}, \code{\link{makeTxDbFromGRanges}},
          and \code{\link{makeTxDbFromGFF}}, for convenient ways to make a
          \link{TxDb} object from UCSC online resources, or from a
          \link[GenomicRanges]{GRanges} object, or from a GFF or GTF file.

    \item The \code{\link[biomaRt]{listMarts}}, \code{\link[biomaRt]{useMart}},
          \code{\link[biomaRt]{listDatasets}}, and
          \code{\link[biomaRt]{listFilters}} functions in the
          \pkg{biomaRt} package.

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
## ---------------------------------------------------------------------
## A. BASIC USAGE
## ---------------------------------------------------------------------

## We can use listDatasets() from the biomaRt package to list the
## datasets available in the "ENSEMBL_MART_ENSEMBL" BioMart database:
library(biomaRt)
listMarts(host="www.ensembl.org")
datasets <- listDatasets(useMart(biomart="ENSEMBL_MART_ENSEMBL",
                                 host="www.ensembl.org"))
head(datasets)
subset(datasets, grepl("elegans", dataset, ignore.case=TRUE))

## Retrieve the full transcript dataset for Worm:
txdb1 <- makeTxDbFromBiomart(dataset="celegans_gene_ensembl")
txdb1

## Retrieve an incomplete transcript dataset for Human:
transcript_ids <- c(
    "ENST00000013894",
    "ENST00000268655",
    "ENST00000313243",
    "ENST00000435657",
    "ENST00000384428",
    "ENST00000478783"
)
txdb2 <- makeTxDbFromBiomart(dataset="hsapiens_gene_ensembl",
                             transcript_ids=transcript_ids)
txdb2  # note that these annotations match the GRCh38 genome assembly

## ---------------------------------------------------------------------
## B. USING A HOST OTHER THAN www.ensembl.org
## ---------------------------------------------------------------------

## A typical use case is to access the "ENSEMBL_MART_ENSEMBL" BioMart
## database on a mirror e.g. on uswest.ensembl.org. A gotcha when
## doing this is that the name of the database on the mirror might
## be different! We can check this with listMarts() from the biomaRt
## package:
listMarts(host="useast.ensembl.org")

## Therefore in addition to setting 'host' to "uswest.ensembl.org" we
## might also need to specify the 'biomart' argument:
txdb3 <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                             dataset="hsapiens_gene_ensembl",
                             transcript_ids=transcript_ids,
                             host="useast.ensembl.org")
txdb3

## ---------------------------------------------------------------------
## C. USING FILTERS
## ---------------------------------------------------------------------

## We can use listFilters() from the biomaRt package to get valid filter
## names:
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl",
                host="www.ensembl.org")
head(listFilters(mart))

## Retrieve transcript dataset for Ensembl gene ENSG00000011198:
my_filter <- list(ensembl_gene_id="ENSG00000011198")
txdb4 <- makeTxDbFromBiomart(dataset="hsapiens_gene_ensembl",
                             filters=my_filter)
txdb4
transcripts(txdb4, columns=c("tx_id", "tx_name", "gene_id"))
transcriptLengths(txdb4)

## ---------------------------------------------------------------------
## D. RETRIEVING CHROMOSOME INFORMATION ONLY
## ---------------------------------------------------------------------

chrominfo <- getChromInfoFromBiomart(dataset="celegans_gene_ensembl")
chrominfo
}
