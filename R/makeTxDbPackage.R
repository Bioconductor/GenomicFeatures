## So to make a package we need a couple things:
## 1) we need a method called makeTxDbPackage (that will take a txdb object)
## 2) we will need a package template to use


## Separate helper function for abbreviating the genus and species name strings
.capitalizeFirstLetter <- function(str){
    paste0(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str)))
}

## And now the name can be any length (as long as there are at least 2 strings)
.abbrevOrganismName <- function(organism){
  spc <- unlist(strsplit(organism, " "))
  if(length(spc)<2){
      stop(strwrap(paste0("Organism should have a genus and species separated",
                          " by a space,")))}
  if(length(spc)==2){
      res <- paste0( substr(spc[[1]], 1, 1), spc[[2]])
  }
  if(length(spc)>2){
      res <- paste0(toupper(substr(spc[[1]], 1, 1)), ## capital of genus
                    spc[[2]],                        ## species
                    ## and any subspecies as camelcase
                    paste0(unlist(lapply(spc[3:length(spc)],
                                         .capitalizeFirstLetter)),
                           collapse=""),
                    collapse="")
  }
  res
}

## simplify DB retrievals from metadata table
.getMetaDataValue <- function(txdb, name){
  con <- dbconn(txdb)
  res <- dbGetQuery(con,
    paste0("SELECT value FROM metadata WHERE name='", name, "'"))[[1]]
  if(!length(res))
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

makePackageName <- function(txdb){
  prefix <- .choosePrefix(txdb)
  con <- dbconn(txdb)
  species <- .abbrevOrganismName(.getMetaDataValue(txdb,'Organism'))
  type <- .getMetaDataValue(txdb,'Data source')
  if(type=="UCSC"){
    genome <- .getMetaDataValue(txdb,'Genome')
    table <- .getMetaDataValue(txdb,'UCSC Table')
    pkgName <- paste(prefix,species,type,genome,table, sep=".")
  }else if(type=="BioMart"){
    db <- .getMetaDataValue(txdb,'BioMart database')
    if(tolower(substr(db, 1, 7)) == "ensembl"){
    dbVer <- .getMetaDataValue(txdb,'BioMart dataset version')
      pkgName <- paste(prefix,species,type,db,dbVer, sep=".")
    }else{
      pkgName <- paste(prefix,species,type,db, sep=".")
    }
  }else{
    type <- gsub(" ",".",type)
    pkgName <- paste(prefix,species,type, sep=".")
  }
  gsub("_","",pkgName)  ## R does not allow underscores in package names
}

.makeObjectName <- function(pkgName){
  strs <- unlist(strsplit(pkgName, "\\."))
  paste(c(strs[2:length(strs)],strs[1]), collapse="_")
}

.getTxDbVersion <- function(txdb){
    type <- .getMetaDataValue(txdb,'Data source')

    if (type=="UCSC") {
        version <- paste(.getMetaDataValue(txdb,'Genome'),
                         "genome based on the",
                         .getMetaDataValue(txdb,'UCSC Table'), "table")
    } else if(type=="BioMart") {
      version <- .getMetaDataValue(txdb,'BioMart database version')
    } else {
      version <- .getMetaDataValue(txdb,'Data source')
    }
    version
}

.normMaintainer <- function(maintainer) {
    maintainer <- as.person(maintainer)
    if (length(maintainer) > 1L) {
        stop("more than one 'maintainer' provided")
    }
    maintainer
}

.getMaintainer <- function(authors) {
    m <- vapply(authors, function(a) "cre" %in% a$role, logical(1L))
    if (sum(m) != 1L) {
        stop("there must be one 'maintainer'")
    }
    maintainer <- authors[m]
    maintainer$role <- list(NULL)
    maintainer$comment <- list(NULL)
    maintainer
}

.mergeMaintainer <- function(authors, maintainer) {
    maintainer <- .normMaintainer(maintainer)
    maintainer$role <- list(union(maintainer$role, "cre"))
    m <- unlist(authors$given) == maintainer$given &
        unlist(authors$family) == maintainer$family
    if (any(m)) {
        authors$role[m] <- list(union(unlist(authors$role[m]), "cre"))
        if (!is.null(maintainer$email)) {
            authors$email[m] <- maintainer$email
        }
    } else {
        authors <- c(authors, maintainer)
    }
    maintainer <- .getMaintainer(authors)
    if (is.null(maintainer$email)) {
        stop("the 'maintainer' must have an email address")
    }
    authors
}

.normAuthor <- function(authors, maintainer) {
    authors <- as.person(authors)
    if (!missing(maintainer)) {
        authors <- .mergeMaintainer(authors, maintainer)
    }
    authors
}



#' Making a TxDb package from annotations available at the UCSC Genome Browser,
#' biomaRt or from another source.
#' 
#' A \link{TxDb} package is an annotation package containing a \link{TxDb}
#' object.
#' 
#' The \code{makeTxDbPackageFromUCSC} function allows the user to make a
#' \link{TxDb} package from transcript annotations available at the UCSC Genome
#' Browser.
#' 
#' The \code{makeTxDbPackageFromBiomart} function allows the user to do the
#' same thing as \code{makeTxDbPackageFromUCSC} except that the annotations
#' originate from biomaRt.
#' 
#' Finally, the \code{makeTxDbPackage} function allows the user to make a
#' \link{TxDb} package directly from a \link{TxDb} object.
#' 
#' \code{makeTxDbPackageFromUCSC} is a convenience function that calls both the
#' \code{\link{makeTxDbFromUCSC}} and the \code{\link{makeTxDbPackage}}
#' functions.  The \code{makeTxDbPackageFromBiomart} follows a similar pattern
#' and calls the \code{\link{makeTxDbFromBiomart}} and
#' \code{\link{makeTxDbPackage}} functions.  \code{supportedMiRBaseBuildValues}
#' is a convenience function that will list all the possible values for the
#' miRBaseBuild argument.  \code{makePackageName} creates a package name from a
#' TxDb object.  This function is also used by OrganismDbi.
#' 
#' @aliases makeTxDbPackage makeTxDbPackageFromUCSC makeFDbPackageFromUCSC
#' makeTxDbPackageFromBiomart supportedMiRBaseBuildValues makePackageName
#' @param version What is the version number for this package?
#' @param maintainer Who is the package maintainer? (must include email to be
#' valid). Should be a \code{\link{person}} object, or something coercible to
#' one, like a string. May be omitted if the \code{author} argument is a
#' \code{person} containing someone with the maintainer role.
#' @param author Who is the creator of this package? Should be a
#' \code{\link{person}} object, or something coercible to one, like a character
#' vector of names. The \code{maintainer} argument will be merged into this
#' list.
#' @param destDir A path where the package source should be assembled.
#' @param license What is the license (and it's version)
#' @param biomart which BioMart database to use.  Get the list of all available
#' BioMart databases with the \code{\link[biomaRt]{listMarts}} function from
#' the biomaRt package. See the details section below for a list of BioMart
#' databases with compatible transcript annotations.
#' @param dataset which dataset from BioMart. For example:
#' \code{"hsapiens_gene_ensembl"}, \code{"mmusculus_gene_ensembl"},
#' \code{"dmelanogaster_gene_ensembl"}, \code{"celegans_gene_ensembl"},
#' \code{"scerevisiae_gene_ensembl"}, etc in the ensembl database.  See the
#' examples section below for how to discover which datasets are available in a
#' given BioMart database.
#' @param genome genome abbreviation used by UCSC and obtained by
#' \code{\link[rtracklayer]{ucscGenomes}()[ , "db"]}.  For example:
#' \code{"hg18"}.
#' @param track name of the UCSC track.  Use
#' \code{supportedUCSCFeatureDbTracks} to get the list of available tracks for
#' a particular genome
#' @param tablename name of the UCSC table containing the transcript
#' annotations to retrieve. Use the \code{\link{supportedUCSCtables}} utility
#' function to get the list of tables known to work with
#' \code{makeTxDbFromUCSC}.
#' @param transcript_ids optionally, only retrieve transcript annotation data
#' for the specified set of transcript ids.  If this is used, then the meta
#' information displayed for the resulting \link{TxDb} object will say 'Full
#' dataset: no'.  Otherwise it will say 'Full dataset: yes'.
#' @param circ_seqs a character vector to list out which chromosomes should be
#' marked as circular.
#' @param filter Additional filters to use in the BioMart query. Must be a
#' named list. An example is \code{filter=as.list(c(source="entrez"))}
#' @param host The host URL of the BioMart. Defaults to
#' https://www.ensembl.org.
#' @param port The port to use in the HTTP communication with the host. This
#' argument has been deprecated. It is handled by \code{useEnsembl} depending
#' on the host input.
#' @param id_prefix Specifies the prefix used in BioMart attributes. For
#' example, some BioMarts may have an attribute specified as
#' \code{"ensembl_transcript_id"} whereas others have the same attribute
#' specified as \code{"transcript_id"}. Defaults to \code{"ensembl_"}.
#' @param columns a named character vector to list out the names and types of
#' the other columns that the downloaded track should have.  Use
#' \code{UCSCFeatureDbTableSchema} to retrieve this information for a
#' particular table.
#' @param url,goldenPath.url use to specify the location of an alternate UCSC
#' Genome Browser.
#' @param chromCol If the schema comes back and the 'chrom' column has been
#' labeled something other than 'chrom', use this argument to indicate what
#' that column has been labeled as so we can properly designate it.  This could
#' happen (for example) with the knownGene track tables, which has no
#' 'chromStart' or 'chromEnd' columns, but which DOES have columns that could
#' reasonably substitute for these columns under particular circumstances.
#' Therefore we allow these three columns to have arguments so that their
#' definition can be re-specified
#' @param chromStartCol Same thing as chromCol, but for renames of 'chromStart'
#' @param chromEndCol Same thing as chromCol, but for renames of 'chromEnd'
#' @param txdb A \link{TxDb} object that represents a handle to a transcript
#' database. This object type is what is returned by \code{makeTxDbFromUCSC},
#' \code{makeTxDbFromUCSC} or \code{makeTxDb}
#' @param taxonomyId By default this value is NA and the organism provided (or
#' inferred) will be used to look up the correct value for this.  But you can
#' use this argument to override that and supply your own valid taxId here
#' @param miRBaseBuild specify the string for the appropriate build Information
#' from mirbase.db to use for microRNAs.  This can be learned by calling
#' \code{supportedMiRBaseBuildValues}.  By default, this value will be set to
#' \code{NA}, which will inactivate the \code{microRNAs} accessor.
#' @param pkgname By default this value is NULL and does not need to be filled
#' in (a package name will be generated for you).  But if you override this
#' value, then the package and it's object will be instead named after this
#' value.  Be aware that the standard rules for package names will apply, (so
#' don't include spaces, underscores or dashes)
#' @param provider If not given, a default is taken from the 'Data source'
#' field of the metadata table.
#' @param providerVersion If not given, a default is taken from one of 'UCSC
#' table', 'BioMart version' or 'Data source' fields of the metadata table.
#' @return A \link{TxDb} object.
#' @author M. Carlson
#' @seealso \code{\link{makeTxDbFromUCSC}}, \code{\link{makeTxDbFromBiomart}},
#' \code{\link{makeTxDb}}, \code{\link[rtracklayer]{ucscGenomes}}
#' @examples
#' 
#' ## First consider relevant helper/discovery functions:
#' ## Get the list of tables known to work with makeTxDbPackageFromUCSC():
#' supportedUCSCtables(genome="hg19")
#' 
#' ## Can also list all the possible values for the miRBaseBuild argument:
#' supportedMiRBaseBuildValues()
#' 
#' ## Next are examples of actually building a package:
#' \dontrun{
#' ## Makes a transcript package for Yeast from the ensGene table at UCSC:
#' makeTxDbPackageFromUCSC(version="0.01",
#'                         maintainer="Some One <so@someplace.org>",
#'                         author="Some One <so@someplace.com>",
#'                         genome="sacCer2",
#'                         tablename="ensGene")
#' 
#' ## Makes a transcript package from Human by using biomaRt and limited to a
#' ## small subset of the transcripts.
#' transcript_ids <- c(
#'     "ENST00000400839",
#'     "ENST00000400840",
#'     "ENST00000478783",
#'     "ENST00000435657",
#'     "ENST00000268655",
#'     "ENST00000313243",
#'     "ENST00000341724")
#' 
#' makeTxDbPackageFromBiomart(version="0.01",
#'                            maintainer="Some One <so@someplace.org>",
#'                            author="Some One <so@someplace.com>",
#'                            transcript_ids=transcript_ids)
#' 
#' }
#' 
#' 
#' @export makeTxDbPackage
makeTxDbPackage <- function(txdb,
                            version,
                            maintainer,
                            author,
                            destDir=".",
                            license="Artistic-2.0",
                            pkgname=NULL,
                            provider=NULL,
                            providerVersion=NULL){
   ## every package has a name We will generate this according to a heuristic
   if (is.null(pkgname))
       pkgname <- makePackageName(txdb)
   if (is.null(provider))
       provider <- .getMetaDataValue(txdb,'Data source')
   else if (!isSingleString(provider))
       stop("'provider' must be a single string")
   if (is.null(providerVersion))
       providerVersion <- .getTxDbVersion(txdb)
   else if (!isSingleString(as.character(providerVersion)))
       stop("'providerVersion' must be a single string")
   dbType <- .getMetaDataValue(txdb,'Db type')
   authors <- .normAuthor(author, maintainer)

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
    AUTHOR=paste(authors, collapse=", "),
    MAINTAINER=as.character(.getMaintainer(authors)),
    GFVERSION=.getMetaDataValue(txdb,
      'GenomicFeatures version at creation time'),
    LIC=license,
    DBTYPE= dbType,
    ORGANISM=.getMetaDataValue(txdb,'Organism'),
    SPECIES=.getMetaDataValue(txdb,'Organism'),
    PROVIDER=provider,
    PROVIDERVERSION=providerVersion,
    RELEASEDATE=.getMetaDataValue(txdb,'Creation time'),
    SOURCEURL=.getMetaDataValue(txdb,'Resource URL'),
    ORGANISMBIOCVIEW=gsub(" ","_",.getMetaDataValue(txdb,'Organism')),
    ## For now: keep conventional object names
    TXDBOBJNAME=pkgname ## .makeObjectName(pkgname)
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
   createPackage(pkgname=pkgname, destinationDir=destDir,
                 originDir=template_path, symbolValues=symvals)
   ## then copy the contents of the database into the extdata dir
   db_path <- file.path(destDir, pkgname, "inst", "extdata",
                        paste(pkgname,"sqlite",sep="."))
   if (!dir.exists(dirname(db_path)))
       dir.create(dirname(db_path), recursive=TRUE)
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
  circ_seqs=NULL,
  url="http://genome.ucsc.edu/cgi-bin/",
  goldenPath.url=getOption("UCSC.goldenPath.url"),
  taxonomyId=NA,
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
    if(!isSingleString(goldenPath.url)){
        stop("'goldenPath.url' must be supplied as a single element",
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
                             goldenPath.url=goldenPath.url,
                             taxonomyId=taxonomyId,
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
  biomart="ENSEMBL_MART_ENSEMBL",
  dataset="hsapiens_gene_ensembl",
  transcript_ids=NULL,   ## optional
  circ_seqs=NULL,
  filter=NULL,
  id_prefix="ensembl_",
  host="https://www.ensembl.org",
  port,
  taxonomyId=NA,
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
    if (!missing(port))
        warning("The 'port' argument is deprecated and will be ignored.")
    ## Make the DB
    txdb <- makeTxDbFromBiomart(biomart=biomart,
                                dataset=dataset,
                                transcript_ids=transcript_ids,
                                circ_seqs=circ_seqs,
                                filter=filter,
                                id_prefix=id_prefix,
                                host=host,
                                taxonomyId=taxonomyId,
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
    goldenPath.url=getOption("UCSC.goldenPath.url"),
    chromCol=NULL,
    chromStartCol=NULL,
    chromEndCol=NULL,
    taxonomyId=NA){
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
    if(!isSingleString(goldenPath.url)){
        stop("'goldenPath.url' must be supplied as a single element",
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
                                 goldenPath.url=goldenPath.url,
                                 chromCol=chromCol,
                                 chromStartCol=chromStartCol,
                                 chromEndCol=chromEndCol,
                                 taxonomyId=taxonomyId)

    ## Make the Package (recycle functions and templates from txdb)
    makeTxDbPackage(fdb,
                    version=version,
                    maintainer=maintainer,
                    author=author,
                    destDir=destDir,
                    license=license)
}

