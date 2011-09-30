###
### Load the txdb object for them when the package is loaded.
###

.onLoad <- function(libname, pkgname)
{
  txdb <- loadDb(system.file("extdata", paste(pkgname,
    ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)
  objname <- "@TXDBOBJNAME@"
  ns <- asNamespace(pkgname)
  assign(objname, txdb, envir=ns)
  namespaceExport(ns, objname)
}

