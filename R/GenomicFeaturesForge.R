###

read_knownGene_table <- function(file)
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
    read.table(file, sep="\t", col.names=names(COL2CLASS),
               colClasses=COL2CLASS, check.names=FALSE)
}

