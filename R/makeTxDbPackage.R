## So to make a package we need a couple things:
## 1) we need a method called makeTxDbPackage (that will take a txdb object)
## 2) we will need a package template to use


## Separate helper function for abbreviating the genus and species name strings
.abbrevOrganismName <- function(organism){
  spc <- unlist(strsplit(organism, " "))
  ## this assumes a binomial nomenclature has been maintained.
  paste0( substr(spc[[1]], 1, 1), spc[[2]])
}

## simplify DB retrievals from metadata table
.getMetaDataValue <- function(txdb, name){
  con <- AnnotationDbi:::dbconn(txdb)
  res <- dbGetQuery(con, 
    paste0("SELECT value FROM metadata WHERE name='", name, "'"))[[1]]  
  if(!is.character(res))
      stop("metadata table missing a value for '", name, "'")
  res
}


## helper functions
.choosePrefix <- function(txdb){
  pkgType <- .getMetaDataValue(txdb,'Db type')
  if(pkgType %in% c("TranscriptDb", "TxDb")){
    prefix <- "TxDb"
  }else if(pkgType == "FeatureDb"){
    prefix <- "FDb"
  }
  prefix
}

.makePackageName <- function(txdb){
  prefix <- .choosePrefix(txdb)
  con <- AnnotationDbi:::dbconn(txdb)
  species <- .abbrevOrganismName(.getMetaDataValue(txdb,'Organism'))  
  type <- .getMetaDataValue(txdb,'Data source')
  if(type=="UCSC"){
    genome <- .getMetaDataValue(txdb,'Genome')
    table <- .getMetaDataValue(txdb,'UCSC Table')
    pkgName <- paste(prefix,species,type,genome,table, sep=".")
  }else if(type=="BioMart"){
    db <- .getMetaDataValue(txdb,'BioMart database')
    if(db == "ensembl"){
    dbVer <- .getMetaDataValue(txdb,'BioMart dataset version')      
      pkgName <- paste(prefix,species,type,db,dbVer, sep=".")      
    }else{
      pkgName <- paste(prefix,species,type,db, sep=".")      
    }
  }else{
    type <- gsub(" ",".",type)
    pkgName <- paste(prefix,species, sep=".")      
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
    version <- paste(.getMetaDataValue(txdb,'Genome'),"genome based on the", 
      .getMetaDataValue(txdb,'UCSC Table'), "table")
  }else if(type=="BioMart"){
    version <- .getMetaDataValue(txdb,'BioMart database version')   
  }else{
    version <-  .getMetaDataValue(txdb,'Data source')   
  }
  version 
}

.getResourceURLInfo <- function(txdb){
    type <- .getMetaDataValue(txdb,'Data source')
    if(type=="UCSC" || type=="BioMart"){
        url <- .getMetaDataValue(txdb,'Resource URL')
    }else{
        url <- .getMetaDataValue(txdb,'Data source')
    }
    url
}


makeTxDbPackage <- function(txdb,
                            version,
                            maintainer,
                            author,
                            destDir=".",
                            license="Artistic-2.0"){
   ## every package has a name We will generate this according to a heuristic
   pkgName <- .makePackageName(txdb)
   dbType <- .getMetaDataValue(txdb,'Db type')
   
   ## there should only be one template
   template_path <- system.file("txdb-template",package="GenomicFeatures")
   ## We need to define some symbols in order to have the 
   ## template filled out correctly.
   symvals <- list(
    PKGTITLE=paste("Annotation package for",dbType,
      "object(s)"),
    PKGDESCRIPTION=paste("Exposes an annotation databases generated from",
      .getMetaDataValue(txdb,'Data source'), "by exposing these as",dbType,
      "objects"),
    PKGVERSION=version,
    AUTHOR=author,
    MAINTAINER=maintainer,
    GFVERSION=.getMetaDataValue(txdb,
      'GenomicFeatures version at creation time'),
    LIC=license,
    DBTYPE= dbType,
    ORGANISM=.getMetaDataValue(txdb,'Organism'),
    SPECIES=.getMetaDataValue(txdb,'Organism'),
    PROVIDER=.getMetaDataValue(txdb,'Data source'),
    PROVIDERVERSION=.getTxDbVersion(txdb),
    RELEASEDATE= .getMetaDataValue(txdb,'Creation time'),
    SOURCEURL= .getResourceURLInfo(txdb),
    ORGANISMBIOCVIEW=gsub(" ","_",.getMetaDataValue(txdb,'Organism')),
    TXDBOBJNAME=pkgName ## .makeObjectName(pkgName)
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
   saveDb(txdb, file=db_path)
}




## Questions: 1) should I set defaults to NULL or use other strategy?
## A: use missing()

## AlSO: 1) put herve quotes around the argument names.
## 2) use IsSingleString (from IRanges) and keep the error IN the
## function (so context is maintained)
## To handle 40 of the same warning: withCallingHandlers()
## elaborate example of how to do a good job of wrangling multiple warnings()

## library(GenomicFeatures);txdb <- makeTxDbPackageFromUCSC(version='1', author="me", maintainer="you',genome = 'mm9', tablename = 'refGene')

## THIS STILL DOESN'T WORK RIGHT!
## IOW THIS FAILS
## txdb <- makeTxDbPackageFromUCSC(version="1.0", genome = 'mm9', tablename = 'refGene')

makeTxDbPackageFromUCSC <- function(
  version,   
  maintainer,
  author,
  destDir=".",
  license="Artistic-2.0",
  genome="hg19",
  tablename="knownGene",
  transcript_ids=NULL,    ## optional
  circ_seqs=DEFAULT_CIRC_SEQS,
  url="http://genome.ucsc.edu/cgi-bin/",
  goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
  miRBaseBuild=NA){
    ## checks
    if(missing(version) || !isSingleString(version)){
        stop("'version' must be supplied as a single element",
             " character vector.")}
    if(missing(version) || !isSingleString(maintainer)){
        stop("'maintainer' must be supplied as a single element",
             " character vector.")}
    if(missing(version) || !isSingleString(author)){
        stop("'author' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(destDir)){
        stop("'destDir' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(license)){
        stop("'license' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(genome)){
        stop("'genome' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(tablename)){
        stop("'tablename' must be supplied as a single element",
             " character vector.")}
    if(!is.character(circ_seqs) || length(circ_seqs)<1){
        stop("'circ_seqs' must be supplied as a named character vector.")}
    if(!isSingleString(url)){
        stop("'url' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(goldenPath_url)){
        stop("'goldenPath_url' must be supplied as a single element",
             " character vector.")}
    if(!isSingleStringOrNA(miRBaseBuild)){
        stop("'miRBaseBuild' must be supplied as a single element",
             " character vector or be NA.")}
    
    ## Make the DB
    txdb <- makeTxDbFromUCSC(genome=genome,
                             tablename=tablename,
                             transcript_ids=transcript_ids,
                             circ_seqs=circ_seqs,
                             url=url,
                             goldenPath_url=goldenPath_url,
                             miRBaseBuild=miRBaseBuild)
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
  transcript_ids=NULL,   ## optional
  circ_seqs=DEFAULT_CIRC_SEQS, 
  filters="",
  id_prefix="ensembl_",
  host="www.biomart.org",
  port=80,
  miRBaseBuild=NA){
    ## checks
    if(missing(version) || !isSingleString(version)){
        stop("'version' must be supplied as a single element",
             " character vector.")}
    if(missing(maintainer) || !isSingleString(maintainer)){
        stop("'maintainer' must be supplied as a single element",
             " character vector.")}
    if(missing(author) || !isSingleString(author)){
        stop("'author' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(destDir)){
        stop("'destDir' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(license)){
        stop("'license' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(biomart)){
        stop("'biomart' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(dataset)){
        stop("'dataset' must be supplied as a single element",
             " character vector.")}
    if(!is.character(circ_seqs) || length(circ_seqs)<1){
        stop("'circ_seqs' must be supplied as a named character vector.")}
    if(!isSingleStringOrNA(miRBaseBuild)){
        stop("'miRBaseBuild' must be supplied as a single element",
             " character vector or be NA.")}
    ## Make the DB
    txdb <- makeTxDbFromBiomart(biomart=biomart,
                                dataset=dataset,
                                transcript_ids=transcript_ids,
                                circ_seqs=circ_seqs,
                                filters=filters,
                                id_prefix=id_prefix,
                                host=host,
                                port=port,
                                miRBaseBuild=miRBaseBuild)
    ## Make the Package
    makeTxDbPackage(txdb,
                    version=version,
                    maintainer=maintainer,
                    author=author,
                    destDir=destDir,
                    license=license)
}



## One for FeatureDB
makeFDbPackageFromUCSC <- function(
    version,
    maintainer,
    author,
    destDir=".",
    license="Artistic-2.0",
    genome="hg19",
    track="tRNAs",
    tablename="tRNAs",
    columns = UCSCFeatureDbTableSchema(genome, track, tablename), 
    url="http://genome.ucsc.edu/cgi-bin/",
    goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
    chromCol=NULL,
    chromStartCol=NULL,
    chromEndCol=NULL){
    ## checks
    if(missing(author) || !isSingleString(version)){
        stop("'version' must be supplied as a single element",
             " character vector.")}
    if(missing(author) || !isSingleString(maintainer)){
        stop("'maintainer' must be supplied as a single element",
             " character vector.")}
    if(missing(author) || !isSingleString(author)){
        stop("'author' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(destDir)){
        stop("'destDir' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(license)){
        stop("'license' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(genome)){
        stop("'genome' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(track)){
        stop("'track' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(tablename)){
        stop("'tablename' must be supplied as a single element",
             " character vector.")}    
    if(!is.character(columns) || length(columns)<1  ||
       length(names(columns))==0 ){
        stop("columns must be supplied as a named character vector.")}
    if(!isSingleString(url)){
        stop("'url' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(goldenPath_url)){
        stop("'goldenPath_url' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(chromCol)){
        stop("'chromCol' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(chromStartCol)){
        stop("'chromStartCol' must be supplied as a single element",
             " character vector.")}
    if(!isSingleString(chromEndCol)){
        stop("'chromEndCol' must be supplied as a single element",
             " character vector.")}
   
    ## make the fdb
    fdb <- makeFeatureDbFromUCSC(genome=genome,
                                 track=track,
                                 tablename=tablename,
                                 columns=columns,
                                 url=url,
                                 goldenPath_url=goldenPath_url,
                                 chromCol=chromCol,
                                 chromStartCol=chromStartCol,
                                 chromEndCol=chromEndCol)
  
    ## Make the Package (recycle functions and templates from txdb)
    makeTxDbPackage(fdb,
                    version=version,
                    maintainer=maintainer,
                    author=author,
                    destDir=destDir,
                    license=license)
}

