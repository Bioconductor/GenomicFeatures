\name{@TXDBOBJNAME@}
\docType{package}

\alias{@PKGNAME@-package}
\alias{@PKGNAME@}
\alias{@TXDBOBJNAME@}


\title{@PKGTITLE@}

\description{
  @PKGDESCRIPTION@
  
  This package loads the above mentioned object.  This TranscriptDb object is an 
  R interface to the prefabricated database contained by this package.
}

\note{
  This Transcript data package was made from resources at @PROVIDER@ on @RELEASEDATE@ and based on the @PROVIDERVERSION@
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
## show the txdb object that is loaded by calling it's name
@PKGNAME@

}

\keyword{package}
\keyword{data}
