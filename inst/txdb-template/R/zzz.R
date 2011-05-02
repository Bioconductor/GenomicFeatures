###
### Load the txdb object for them when the package is loaded.
###

.onLoad <- function(libname, pkgname)
{
  require("GenomicFeatures", quietly=TRUE)
  @TXDBOBJNAME@ <- loadFeatures(system.file("extdata", paste(pkgname,
    ".sqlite",sep=""), package=pkgname, lib.loc=libname))
}

