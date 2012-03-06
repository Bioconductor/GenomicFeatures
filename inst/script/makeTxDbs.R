library(GenomicFeatures)

version = "2.6.4"
maintainer = "Biocore Data Team <bioconductor@r-project.org>"
author = "Marc Carlson"


## Human HG18
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=outDir,
                        genome="hg18",
                        tablename="knownGene")

## Human HG19
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=outDir,
                        genome="hg19",
                        tablename="knownGene",
                        miRBaseBuild="GRCh37")

## Mouse mm9
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=outDir,
                        genome="mm9",
                        tablename="knownGene")



## c elegans enseble trac
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=outDir,
                        genome="ce6",
                        tablename="ensGene")
                        
## Drosophila ensembl track
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=outDir,
                        genome="dm3",
                        tablename="ensGene")
                        
## Arabidopsis biomaRt
makeTxDbPackageFromBiomart(version=version,
                          maintainer=maintainer,
                          author=author,
		 	  destDir=outDir,
                          biomart="plants_mart_12",
                          dataset="athaliana_eg_gene")
                                         

## Rat rn4 ensembl track
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=outDir,
                        genome="rn4",
                        tablename="ensGene")



## yeast ensembl genes is no longer available?
## Yeast sacCer2 ensembl track
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=outDir,
                        genome="sacCer2",
                        tablename="sgdGene")

## makeTxDbPackageFromUCSC(version=version,
##                         maintainer=maintainer,
##                         author=author,
##			   destDir=outDir,
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
## 			      destDir=outDir,
##                            biomart="ensembl",
##                            dataset="hsapiens_gene_ensembl",
##                            transcript_ids=transcript_ids,
##                            miRBaseBuild="GRCh37")
                                         


