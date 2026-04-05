.onUnload <- function(libpath)
{
    library.dynam.unload("GenomicFeatures", libpath)
}

.test <- function() BiocGenerics:::testPackage("GenomicFeatures")

