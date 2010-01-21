### Should we rename the first arg. 'txdb' or just 'x'?
### Also, should we change the signature to dispatch on the first arg. only?

## setGeneric("transcriptsByRanges", signature="ranges",
##            function(txdb, ranges=NULL, restrict="any",
##                     columns=c("tx_id", "tx_name"))
##            standardGeneric("transcriptsByRanges"))

## setGeneric("bindTranscripts", signature="ranges",
##            function(txdb, ranges=NULL, restrict="any",
##                     columns=c("tx_id", "tx_name"))
##            standardGeneric("bindTranscripts"))

## setGeneric("transcripts", signature="vals",
##            function(txdb, vals=NULL, columns=c("tx_id", "tx_name"))
##            standardGeneric("transcripts"))





## setGeneric("exonsByRanges", signature="ranges",
##            function(txdb, ranges=NULL, restrict="any",
##                     columns=c("exon_id", "exon_name"))
##            standardGeneric("exonsByRanges"))

## setGeneric("bindExons", signature="ranges",
##            function(txdb, ranges=NULL, restrict="any",
##                     columns=c("exon_id", "exon_name"))
##            standardGeneric("bindExons"))

## setGeneric("exons", signature="vals",
##            function(txdb, vals=NULL, columns=c("tx_id", "tx_name"))
##            standardGeneric("exons"))
