\name{makeTxDbFromGFF}

\alias{makeTxDbFromGFF}

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
                format=c("auto", "gff3", "gtf"),
                dataSource=NA,
                organism=NA,
                taxonomyId=NA,
                circ_seqs=DEFAULT_CIRC_SEQS,
                chrominfo=NULL,
                miRBaseBuild=NA,
                metadata=NULL,
                dbxrefTag)
}
\arguments{
  \item{file}{
    Input GFF3 or GTF file. Can be a path to a file, or an URL,
    or a connection object, or a \link[rtracklayer]{GFF3File}
    or \link[rtracklayer]{GTFFile} object.
  }
  \item{format}{
    Format of the input file. Accepted values are: \code{"auto"} (the
    default) for auto-detection of the format, \code{"gff3"}, or \code{"gtf"}.
    Use \code{"gff3"} or \code{"gtf"} only if auto-detection failed.
  }
  \item{dataSource}{
    A single string describing the origin of the data file. Please be as
    specific as possible.
  }
  \item{organism}{
    What is the Genus and species of this organism. Please use proper
    scientific nomenclature for example: "Homo sapiens" or "Canis familiaris"
    and not "human" or "my fuzzy buddy". If properly written, this information
    may be used by the software to help you out later.
  }
  \item{taxonomyId}{
    By default this value is NA and the organism provided will be used to
    look up the correct value for this. But you can use this argument to
    override that and supply your own taxonomy id here (which will be
    separately validated). Since providing a valid taxonomy id will not
    require us to look up one based on your organism: this is one way that
    you can loosen the restrictions about what is and isn't a valid value
    for the organism.
  }  
  \item{circ_seqs}{
    A character vector to list out which chromosomes should be marked as
    circular.
  }
  \item{chrominfo}{
    Data frame containing information about the chromosomes. Will be
    passed to the internal call to \code{\link{makeTxDb}}.
    See \code{?\link{makeTxDb}} for more information.
    Alternatively, can be a \link[GenomeInfoDb]{Seqinfo} object.
  }
  \item{miRBaseBuild}{
    Specify the string for the appropriate build Information from mirbase.db
    to use for microRNAs. This can be learned by calling
    \code{supportedMiRBaseBuildValues}. By default, this value will be set to
    \code{NA}, which will inactivate the \code{microRNAs} accessor.
  }
  \item{metadata}{
    A 2-column data frame containing meta information to be included in the
    \link{TxDb} object. See \code{?\link{makeTxDb}} for more information about 
    the format of \code{metadata}.
  }
  \item{dbxrefTag}{
    If not missing, the values in the \code{Dbxref} attribute with the
    specified tag (like \dQuote{GeneID}) are used for the feature names.
  }
}

\details{
  \code{makeTxDbFromGFF} is a convenience function that feeds
  data from the parsed file to the \code{\link{makeTxDbFromGRanges}}
  function.
}

\value{A \link{TxDb} object.}

\author{
  M. Carlson and H. Pagès
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
            dataSource="partial gtf file for Tomatoes for testing",
            organism="Solanum lycopersicum")

## TESTING GTF, this time specifying the chrominfo
gtfFile <- system.file("extdata","GTF_files","Aedes_aegypti.partial.gtf",
                       package="GenomicFeatures")
chrominfo <- data.frame(chrom = c('supercont1.1','supercont1.2'),
                        length=c(5220442, 5300000),
                        is_circular=c(FALSE, FALSE))
metadata <- data.frame(name="Resource URL",
                       value=paste0("ftp://ftp.ensemblgenomes.org/pub/metazoa/",
                                    "release-13/gtf/aedes_aegypti/"))
txdb2 <- makeTxDbFromGFF(file=gtfFile,
             chrominfo=chrominfo,
             dataSource="ensemblgenomes",
             organism="Aedes aegypti",
             metadata=metadata)
}
