### =========================================================================
### saveFeatures(), loadFeatures()
### -------------------------------------------------------------------------


setGeneric("saveFeatures", signature="x",
           function(x, file) standardGeneric("saveFeatures"))

setMethod("saveFeatures", "TranscriptDb",
          function(x, file)
          {	      
	    .Defunct(new="saveDb")
            ## saveDb(x, file)
          }
)

setMethod("saveFeatures", "FeatureDb",
          function(x, file)
          {
            .Defunct(new="saveDb")
            ## saveDb(x, file)
          }
)

## Old loadFeatures will be deprecated, but for now lets just not make any
## current users too unhappy.
loadFeatures <- function(file)
{
    .Defunct(new="loadDb")
    ## loadDb(file)
}

