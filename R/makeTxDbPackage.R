### Everything in this file has moved to txdbmaker!

supportedMiRBaseBuildValues <- function()
{
    call_fun_in_txdbmaker("supportedMiRBaseBuildValues")
}

makePackageName <- function(...)
{
    call_fun_in_txdbmaker("makePackageName", ...)
}

makeTxDbPackage <- function(...)
{
    call_fun_in_txdbmaker("makeTxDbPackage", ...)
}

makeTxDbPackageFromUCSC <- function(...)
{
    args <- list(...)
    if ("url" %in% names(args)) {
        .Deprecated(msg="'url' argument is deprecated and was ignored")
        args$url <- NULL
    }
    do.call(call_fun_in_txdbmaker, c(list("makeTxDbPackageFromUCSC"), args))
}

makeFDbPackageFromUCSC <- function(...)
{
    call_fun_in_txdbmaker("makeFDbPackageFromUCSC", ...)
}

makeTxDbPackageFromBiomart <- function(...)
{
    call_fun_in_txdbmaker("makeTxDbPackageFromBiomart", ...)
}

