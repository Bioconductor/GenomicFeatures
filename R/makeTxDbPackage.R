## So to make a package we need a couple things:
## 1) we need a method called makeTxDbPackage (that will take a txdb object)
## 2) we will need a package template to use


## Separate helper function for abbreviating the genus and species name strings
.abbrevSpeciesName <- function(speciesName){
  spc <- unlist(strsplit(speciesName, " "))
  ## this assumes a binomial nomenclature has been maintained.
  paste( substr(spc[[1]], 1, 1), spc[[2]], sep="")
}



## helper functions
.makePackageName <- function(txdb){
  con <- txdbConn(txdb)
  species <- .abbrevSpeciesName(dbGetQuery(con, 
    "SELECT value FROM metadata WHERE name='Genus and Species'")[[1]])  
  type <- dbGetQuery(con, 
    "SELECT value FROM metadata WHERE name='Data source'")[[1]]
  if(type=="UCSC"){
    genome <- dbGetQuery(con, 
    "SELECT value FROM metadata WHERE name='Genome'")[[1]]
    table <- dbGetQuery(con, 
    "SELECT value FROM metadata WHERE name='UCSC Table'")[[1]]
    pkgName <- paste("TxDb",species,type,genome,table, sep=".")
  }else if(type=="BioMart"){
    db <- dbGetQuery(con, 
    "SELECT value FROM metadata WHERE name='BioMart database'")[[1]]
    if(db == "ensembl"){
    dbVer <- dbGetQuery(con, 
    "SELECT value FROM metadata WHERE name='BioMart dataset version'")[[1]]      
      pkgName <- paste("TxDb",species,type,db,dbVer, sep=".")      
    }else{
      pkgName <- paste("TxDb",species,type,db, sep=".")      
    }
  }
  pkgName  
}

.makeTxDbPackage <- function(txdb,
		                         destDir=".",
			                       maintainer,
			                       version){
   ## every package has a name We will generate this according to a heuristic
      pkgName <- .makePackageName(txdb)
   ## 

   createPackage(pkgname=pkgName,
		 destinationDir=destDir,
                 originDir=template_path,
                 symbolValues=symvals)

}
