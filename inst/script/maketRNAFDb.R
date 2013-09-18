## It all has to start with an FDB package
library(GenomicFeatures)
makeFDbPackageFromUCSC(version="1.0.0",
                 maintainer="Bioconductor Package Maintainer <maintainer@bioconductor.org>",
                                author="Marc Carlson",
                                genome="hg19",
                                track="tRNAs",
                                tablename="tRNAs")

## Now I must make more FBDs and save their SQLite files to this package
path <- "FDb.Hsapiens.UCSC.hg19.tRNAs/inst/extdata/"

species <- c(hg18 = "Hsapiens",
             mm9  = "Mmusculus",
             mm10  = "Mmusculus",
             rn4  = "Rnorvegicus",
             rn5  = "Rnorvegicus")

for(i in seq_len(length(species))){
  fdb <- makeFeatureDbFromUCSC(genome=names(species[i]),
                               track="tRNAs",
                               tablename="tRNAs")
  name <- paste("FDb.",species[i],".UCSC.",names(species[i]),
                ".tRNAs.sqlite",sep="")
  saveDb(fdb, file=paste(path,name,sep=""))
}


################################################################################
## Then change what we call it:
file.rename("FDb.Hsapiens.UCSC.hg19.tRNAs","FDb.UCSC.tRNAs")
## Then rename it in DESCRIPTION
fl <- file("FDb.UCSC.tRNAs/DESCRIPTION")
lns <- readLines(fl)
lns <- gsub("FDb.Hsapiens.UCSC.hg19.tRNAs","FDb.UCSC.tRNAs",lns)
writeLines(lns, con = fl)
## Make some edits to the manual page
fl <- file("FDb.UCSC.tRNAs/man/package.Rd")
lns <- readLines(fl)
lns <- gsub("FDb.Hsapiens.UCSC.hg19.tRNAs","FDb.UCSC.tRNAs",lns)
## except for the last one which should always be on line 41
lns[[41]] <-  sub("FDb.UCSC.tRNAs","FDb.Hsapiens.UCSC.hg19.tRNAs",lns[[41]])
lns <- c(lns[1:8],
         "\\alias{FDb.Hsapiens.UCSC.hg19.tRNAs}",
         "\\alias{FDb.Hsapiens.UCSC.hg18.tRNAs}",
         "\\alias{FDb.Mmusculus.UCSC.mm9.tRNAs}",
         "\\alias{FDb.Mmusculus.UCSC.mm10.tRNAs}",
         "\\alias{FDb.Rnorvegicus.UCSC.rn4.tRNAs}",
         "\\alias{FDb.Rnorvegicus.UCSC.rn5.tRNAs}",
         lns[9:length(lns)])
writeLines(lns, con = fl)







## OLDE Stuff follows:

## FDb..UCSC..tRNAs <- makeFeatureDbFromUCSC(genome="",
##                              track="tRNAs",
##                              tablename="tRNAs")
## save(, file=paste(path,".sqlite",sep=""))


## FDb.Hsapiens.UCSC.hg18.tRNAs <- makeFeatureDbFromUCSC(genome="hg18",
##                              track="tRNAs",
##                              tablename="tRNAs")
## saveDb(FDb.Hsapiens.UCSC.hg18.tRNAs,
##      file=paste(path,"FDb.Hsapiens.UCSC.hg18.tRNAs.sqlite",sep=""))


## FDb.Mmusculus.UCSC.mm9.tRNAs <- makeFeatureDbFromUCSC(genome="mm9",
##                              track="tRNAs",
##                              tablename="tRNAs")
## saveDb(FDb.Mmusculus.UCSC.mm9.tRNAs,
##      file=paste(path,"FDb.Mmusculus.UCSC.mm9.tRNAs.sqlite",sep=""))


## FDb.Rnorvegicus.UCSC.mm9.tRNAs <- makeFeatureDbFromUCSC(genome="rn4",
##                              track="tRNAs",
##                              tablename="tRNAs")
## saveDb(FDb.Rnorvegicus.UCSC.mm9.tRNAs,
##      file=paste(path,"FDb.Rnorvegicus.UCSC.mm9.tRNAs.sqlite",sep=""))

################################################################################
## And then I have to totally change the way the /R/zzz.R looks - vectorized
## for all.

## Should look more like this
## .onLoad <- function(libname, pkgname)
## {
##   ns <- asNamespace(pkgname)
##   path <- system.file("extdata", package=pkgname, lib.loc=libname)
##   files <- dir(path)
##   for(i in seq_len(length(files))){
##     fdb <- loadDb(system.file("extdata", files[[i]], package=pkgname, 
##                   lib.loc=libname),packageName=pkgname)
##     objname <- sub(".sqlite$","",files[[i]])
##     assign(objname, fdb, envir=ns)
##     namespaceExport(ns, objname)
##   }
## }



