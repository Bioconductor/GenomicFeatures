\name{@TXDBOBJNAME@}
\docType{package}

\alias{@PKGNAME@-package}
\alias{@PKGNAME@}
\alias{@TXDBOBJNAME@}


\title{@PKGTITLE@}

\description{
  @PKGDESCRIPTION@
}

\details{
  TODO Details about the package???
}

\note{
  This Transcript data package was made from resources at @PROVIDER@ on @RELEASEDATE@ and based on the @PROVIDERVERSION@
}

\author{@AUTHOR@}

\seealso{

}

\examples{

if ("AGAPS" \%in\% masknames(@TXDBOBJNAME@)) {

#  ## Check that the assembly gaps contain only Ns:
#  checkOnlyNsInGaps <- function(seq)
#  {
#    ## Replace all masks by the inverted AGAPS mask
#    masks(seq) <- gaps(masks(seq)["AGAPS"])
#    unique_letters <- uniqueLetters(seq)
#    if (any(unique_letters != "N"))
#        stop("assembly gaps contain more than just Ns")
#  }

#  ## A message will be printed each time a sequence is removed
#  ## from the cache:
#  options(verbose=TRUE)

#  for (seqname in seqnames(@TXDBOBJNAME@)) {
#    cat("Checking sequence", seqname, "... ")
#    seq <- TXDBOBJNAME@[[seqname]]
#    checkOnlyNsInGaps(seq)
#    cat("OK\n")
#  }
}

### See the GenomeSearching vignette in the BSgenome software
### package for some examples of genome-wide motif searching using
### Biostrings and the BSgenome data packages:
#if (interactive())
#    vignette("GenomeSearching", package="BSgenome")
#}

\keyword{package}
\keyword{data}
