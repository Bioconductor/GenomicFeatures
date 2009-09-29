###
### See http://genome.ucsc.edu/FAQ/FAQtracks#tracks1 for why the start
### coordinates found in the database files from UCSC need to be translated
### from zero-based to one-based.
###

### Can be used with 'shift=0L' just to trim the trailing comma found
### in the UCSC multivalued fields like 'exonStarts' or 'exonEnds'.
.shift.coordsInMultivaluedField <- function(x, shift)
{
    if (length(x) == 0L)
        return(character(0))
    sapply(strsplit(x, ",", fixed=TRUE),
           function(starts)
               paste(as.character(as.integer(starts) + shift), collapse=",")
    )
}

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
    return(ans)
}

