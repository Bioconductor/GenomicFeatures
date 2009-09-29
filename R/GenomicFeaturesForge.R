###
### See http://genome.ucsc.edu/FAQ/FAQtracks#tracks1 for why the start
### coordinates found in the database files from UCSC need to be translated
### from zero-based to one-based.
###

### Can be used with 'shift=0L' just to trim the trailing comma found
### in the UCSC multivalued fields like 'exonStarts' or 'exonEnds'.
.shift.coordsInMultivaluedField <- function(x, shift)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    if (length(x) == 0L)
        return(character(0))
    sapply(strsplit(x, ",", fixed=TRUE),
           function(starts)
               paste(as.character(as.integer(starts) + shift), collapse=",")
    )
}

.replaceEmptyStringWithNA <- function(x)
{
    if (!is.character(x))
        stop("'x' must be a character vector")
    if (length(x) == 0L)
        return(character(0))
    x[x == ""] <- NA_character_
    return(x)
}

### Read UCSC knownGene.txt or knownGeneOld3.txt database files.
read_knownGene_table <- function(file, translate.starts=TRUE)
{
    COL2CLASS <- c(
        `name`="character",
        `chrom`="character",
        `strand`="character",
        `txStart`="integer",
        `txEnd`="integer",
        `cdsStart`="integer",
        `cdsEnd`="integer",
        `exonCount`="integer",
        `exonStarts`="character",
        `exonEnds`="character",
        `proteinID`="character",  # use "NULL" to drop
        `alignID`="character"     # use "NULL" to drop
    )
    ans <- read.table(file, sep="\t", col.names=names(COL2CLASS),
                      colClasses=COL2CLASS, check.names=FALSE)
    if (translate.starts) {
        ans$txStart <- ans$txStart + 1L
        ans$cdsStart <- ans$cdsStart + 1L
        ans$exonStarts <- .shift.coordsInMultivaluedField(ans$exonStarts, 1L)
    } else {
        ans$exonStarts <- .shift.coordsInMultivaluedField(ans$exonStarts, 0L)
    }
    ans$exonEnds <- .shift.coordsInMultivaluedField(ans$exonEnds, 0L)
    ans$proteinID <- .replaceEmptyStringWithNA(ans$proteinID)
    return(ans)
}

### Read UCSC cpgIslandExt.txt database file.
read_cpgIslandExt_table <- function(file, translate.starts=TRUE)
{
    COL2CLASS <- c(
        `bin`="NULL",              # dropped
        `chromosome`="character",  # UCSC field is 'chrom'
        `start`="integer",         # UCSC field is 'chromStart'
        `end`="integer",           # UCSC field is 'chromEnd'
        `ID`="character",          # UCSC field is 'name'
        `length`="NULL",           # dropped
        `cpgNum`="NULL",           # dropped
        `gcNum`="NULL",            # dropped
        `perCpg`="NULL",           # dropped
        `perGc`="NULL",            # dropped
        `obsExp`="NULL"            # dropped
    )
    ans <- read.table(file, sep="\t", col.names=names(COL2CLASS),
                      colClasses=COL2CLASS, check.names=FALSE)
    if (translate.starts) {
        ans$start <- ans$start + 1L
    }
    return(ans)
}

