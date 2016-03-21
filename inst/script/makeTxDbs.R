### This script is called from makeTerminalDBPkgs.R
### which is part of the release build code base. The
### version is set in the user workspace/env and not 
### explicitly passed.

library(GenomicFeatures)

maintainer = "Bioconductor Package Maintainer <maintainer@bioconductor.org>"
author = "Marc Carlson"

## Human HG18
cat("building hg18 \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="hg18",
                        tablename="knownGene")

## Human HG19
cat("building hg19 \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="hg19",
                        tablename="knownGene",
                        miRBaseBuild="GRCh37")

## Human HG19 lincRNAs
cat("building hg19 lincRNAs \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="hg19",
                        tablename="lincRNAsTranscripts",
                        miRBaseBuild="GRCh37")


## Looks like UCSC took away their 'knownGene' track for hg38...
## ## Human HG38
## makeTxDbPackageFromUCSC(version=version,
##                         maintainer=maintainer,
##                         author=author,
## 			destDir=TxDbOutDir,
##                         genome="hg38",
##                         tablename="knownGene")



## Mouse mm9
cat("building mm9 \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="mm9",
                        tablename="knownGene")

## Mouse mm10
cat("building mm10 \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="mm10",
                        tablename="ensGene")

## Mouse mm10 knownGene
cat("building mm10 knownGene \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="mm10",
                        tablename="knownGene")


## c elegans enseble track
cat("building ce6 \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="ce6",
                        tablename="ensGene")
                        
## Drosophila ensembl track
cat("building dm3 \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="dm3",
                        tablename="ensGene")


## There is apparently a problem with this resource here:
## It looks like the for the splicing is has missing fields...
## Arabidopsis biomaRt
cat("building plants_mart_28 \n")
makeTxDbPackageFromBiomart(version=version,
                          maintainer=maintainer,
                          author=author,
		 	  destDir=TxDbOutDir,
                          biomart="plants_mart_28",
                          dataset="athaliana_eg_gene")
                                         

## Rat rn4 ensembl track
cat("building rn4 \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="rn4",
                        tablename="ensGene")

## Rat rn5 ensembl track
cat("building rn5 \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="rn5",
                        tablename="refGene")



## yeast ensembl genes is no longer available?
## Yeast sacCer2 ensembl track
cat("building sacCer2 \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="sacCer2",
                        tablename="sgdGene")

cat("building sacCer3 \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="sacCer3",
                        tablename="sgdGene")

## makeTxDbPackageFromUCSC(version=version,
##                         maintainer=maintainer,
##                         author=author,
##			   destDir=TxDbOutDir,
##                         genome="sacCer2",
##                         tablename="ensGene")




## ## ## Human biomaRt (just used for testing)
## transcript_ids <- c(
##                     "ENST00000268655",
##                     "ENST00000313243",
##                     "ENST00000341724",
##                     "ENST00000400839",
##                     "ENST00000400840",
##                     "ENST00000435657",
##                     "ENST00000478783"
##                     )
## ## #debug(GenomicFeatures:::.prepareBiomartMetadata)
## makeTxDbPackageFromBiomart(version=version,
##                            maintainer=maintainer,
##                            author=author,
## 			      destDir=TxDbOutDir,
##                            biomart="ensembl",
##                            dataset="hsapiens_gene_ensembl",
##                            transcript_ids=transcript_ids,
##                            miRBaseBuild="GRCh37")
                                         


