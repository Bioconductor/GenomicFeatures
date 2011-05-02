###
### Load the txdb object for them when the package is loaded.
###

.onLoad <- function(libname, pkgname)
{
  txdb <- loadFeatures(system.file("extdata", paste(pkgname,
    ".sqlite",sep=""), package=pkgname, lib.loc=libname))
  objname <- "@TXDBOBJNAME@"
  ns <- asNamespace(pkgname)
  assign(objname, txdb, envir=ns)
  namespaceExport(ns, objname)
}

