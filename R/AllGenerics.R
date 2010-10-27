##TODO: move the accessor generics into here!

setGeneric("saveFeatures", signature="x",
           function(x, file) standardGeneric("saveFeatures"))


## TODO: loadFeatures needs an argument that it can dispatch on.  As things
## stand now, it doesn't know which method to call.  I should probably demote
## loadFeatures to be a function again, and just have it look inside of the
## .sqlite file and then decide what kind of object to make...
## setGeneric("loadFeatures", signature="file",
##            function(file) standardGeneric("loadFeatures"))
