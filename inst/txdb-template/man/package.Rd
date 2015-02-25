\name{@TXDBOBJNAME@}
\docType{package}

\alias{@PKGNAME@-package}
\alias{@PKGNAME@}
\alias{@TXDBOBJNAME@}


\title{@PKGTITLE@}

\description{  
  This package loads one or more @DBTYPE@ objects.  Such @DBTYPE@
  objects are an R interface to prefabricated databases contained by
  this package.

  The names of any objects exposed by this package indicate the origin and
  resources exposed.  So for example TxDb.Hsapiens.UCSC.hg19.knownGene
  would be a TxDb object, of Homo sapiens data from UCSC build
  hg19 based on the knownGene Track.
}

\note{
  This data package was made from resources at @PROVIDER@ on
  @RELEASEDATE@ and based on the @PROVIDERVERSION@
}

\author{@AUTHOR@}


\seealso{
  \link[GenomicFeatures]{transcripts},
  \link[GenomicFeatures]{transcriptsBy}
}

\examples{
## load the library
library(@PKGNAME@)
## list the contents that are loaded into memory
ls('package:@PKGNAME@')
## show the db object that is loaded by calling it's name
@PKGNAME@

}

\keyword{package}
\keyword{data}
