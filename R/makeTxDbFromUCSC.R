### =========================================================================
### makeTxDbFromUCSC()
### -------------------------------------------------------------------------

### Everything in this file has moved to txdbmaker!

supportedUCSCtables <- function(...)
{
    args <- list(...)
    if ("url" %in% names(args)) {
        .Deprecated(msg="'url' argument is deprecated and was ignored")
        args$url <- NULL
    }
    do.call(call_fun_in_txdbmaker, c(list("supportedUCSCtables"), args))
}

browseUCSCtrack <- function(...)
{
    call_fun_in_txdbmaker("browseUCSCtrack", ...)
}

makeTxDbFromUCSC <- function(...)
{
    args <- list(...)
    if ("url" %in% names(args)) {
        .Deprecated(msg="'url' argument is deprecated and was ignored")
        args$url <- NULL
    }
    do.call(call_fun_in_txdbmaker, c(list("makeTxDbFromUCSC"), args))
}

