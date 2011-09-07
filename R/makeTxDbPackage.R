## So to make a package we need a couple things:
## 1) we need a method called makeTxDbPackage (that will take a txdb object)
## 2) we will need a package template to use


## Separate helper function for abbreviating the genus and species name strings
.abbrevSpeciesName <- function(speciesName){
  spc <- unlist(strsplit(speciesName, " "))
  ## this assumes a binomial nomenclature has been maintained.
  paste( substr(spc[[1]], 1, 1), spc[[2]], sep="")
}

## simplify DB retrievals from metadata table
.getMetaDataValue <- function(txdb, name){
  con <- AnnotationDbi:::dbConn(txdb)
  res <- dbGetQuery(con, 
    paste("SELECT value FROM metadata WHERE name='",
      name,"'", sep=""))[[1]]  
  if(!is.character(res))error("Your metadata table is missing a value for:",
    name,".")
  res
}


## helper functions
.makePackageName <- function(txdb){
  con <- AnnotationDbi:::dbConn(txdb)
  species <- .abbrevSpeciesName(.getMetaDataValue(txdb,'Genus and Species'))  
  type <- .getMetaDataValue(txdb,'Data source')
  if(type=="UCSC"){
    genome <- .getMetaDataValue(txdb,'Genome')
    table <- .getMetaDataValue(txdb,'UCSC Table')
    pkgName <- paste("TxDb",species,type,genome,table, sep=".")
  }else if(type=="BioMart"){
    db <- .getMetaDataValue(txdb,'BioMart database')
    if(db == "ensembl"){
    dbVer <- .getMetaDataValue(txdb,'BioMart dataset version')      
      pkgName <- paste("TxDb",species,type,db,dbVer, sep=".")      
    }else{
      pkgName <- paste("TxDb",species,type,db, sep=".")      
    }
  }
  gsub("_","",pkgName)  ## R does not allow underscores in package names
}

.makeObjectName <- function(pkgName){
  strs <- unlist(strsplit(pkgName, "\\."))
  paste(c(strs[2:length(strs)],strs[1]), collapse="_")
}


.getTxDbVersion <- function(txdb){
  type <- .getMetaDataValue(txdb,'Data source')
  if(type=="UCSC"){  
    version <- paste(.getMetaDataValue(txdb,'Genome'),"genome base on the", 
      .getMetaDataValue(txdb,'UCSC Table'), "table")
  }else if(type=="BioMart"){
    version <- .getMetaDataValue(txdb,'BioMart database version')   
  }
  version 
}


makeTxDbPackage <- function(txdb,
                            version,
			                      maintainer,
                            author,
  	                        destDir=".",
                            license="Artistic-2.0"){
   ## every package has a name We will generate this according to a heuristic
   pkgName <- .makePackageName(txdb)

   ## there should only be one template
   template_path <- system.file("txdb-template",package="GenomicFeatures")
   ## We need to define some symbols in order to have the 
   ## template filled out correctly.
   symvals <- list(
    PKGTITLE=paste("Annotation package for the",.makeObjectName(pkgName),
      "object"),
    PKGDESCRIPTION=paste("Contains the",.makeObjectName(pkgName),"object",
      "annotation database as generated from",
      .getMetaDataValue(txdb,'Data source')),
    PKGVERSION=version,
    AUTHOR=author,
    MAINTAINER=maintainer,
    GFVERSION=.getMetaDataValue(txdb,
      'GenomicFeatures version at creation time'),
    LIC=license,
    ORGANISM=.getMetaDataValue(txdb,'Genus and Species'),
    SPECIES=.getMetaDataValue(txdb,'Genus and Species'),
    PROVIDER=.getMetaDataValue(txdb,'Data source'),
    PROVIDERVERSION=.getTxDbVersion(txdb),
    RELEASEDATE= .getMetaDataValue(txdb,'Creation time'),
    SOURCEURL= .getMetaDataValue(txdb,'Resource URL'),
    ORGANISMBIOCVIEW=gsub(" ","_",.getMetaDataValue(txdb,'Genus and Species')),
    TXDBOBJNAME=.makeObjectName(pkgName)
   )
   ## Should never happen
   if (any(duplicated(names(symvals)))) {
       str(symvals)
       stop("'symvals' contains duplicated symbols")
   }
   ## All symvals should by single strings (non-NA)
   is_OK <- sapply(symvals, isSingleString)
   if (!all(is_OK)) {
       bad_syms <- paste(names(is_OK)[!is_OK], collapse=", ")
       stop("values for symbols ", bad_syms, " are not single strings")
   }
   createPackage(pkgname=pkgName,
		             destinationDir=destDir,
                 originDir=template_path,
                 symbolValues=symvals)
   ## then copy the contents of the database into the extdata dir
   db_path <- file.path(destDir, pkgName, "inst", "extdata", 
     paste(pkgName,"sqlite",sep="."))
   saveFeatures(txdb, file=db_path)
}


## wrapper functions to make BOTH the transcriptDb AND also package it up
## One for UCSC
makeTxDbPackageFromUCSC <- function(
  version,
  maintainer,
  author,
  destDir=".",
  license="Artistic-2.0",
  genome="hg19",
  tablename="knownGene",
  transcript_ids=NULL,
  circ_seqs=DEFAULT_CIRC_SEQS,
  url="http://genome.ucsc.edu/cgi-bin/",
  goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath"){
    ## Make the DB
    txdb <- makeTranscriptDbFromUCSC(genome=genome,
                                     tablename=tablename,
                                     transcript_ids=transcript_ids,
                                     circ_seqs=circ_seqs,
                                     url=url,
                                     goldenPath_url=goldenPath_url)
    ## Make the Package
    makeTxDbPackage(txdb,
                    version=version,
                    maintainer=maintainer,
                    author=author,
                    destDir=destDir,
                    license=license)
}

## One for biomaRt
makeTxDbPackageFromBiomart <- function(
  version,
  maintainer,
  author,
  destDir=".",
  license="Artistic-2.0",
  biomart="ensembl",
  dataset="hsapiens_gene_ensembl",
  transcript_ids=NULL,
  circ_seqs=DEFAULT_CIRC_SEQS){
    ## Make the DB
    txdb <- makeTranscriptDbFromBiomart(biomart=biomart,
                                        dataset=dataset,
                                        transcript_ids=transcript_ids,
                                        circ_seqs=circ_seqs)
    ## Make the Package
    makeTxDbPackage(txdb,
                    version=version,
                    maintainer=maintainer,
                    author=author,
                    destDir=destDir,
                    license=license)
}

