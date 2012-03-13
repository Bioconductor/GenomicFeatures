library(GenomicFeatures)

## Best if I pass this down (from makeTerminal)
## version = "2.7.0" ## best if I pass this down (from makeTerminal)
## TxDbOutDir = paste(outDir,"_TxDbs",sep="")  

maintainer = "Biocore Data Team <bioconductor@r-project.org>"
author = "Marc Carlson"

## Human HG18
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="hg18",
                        tablename="knownGene")

## Human HG19
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="hg19",
                        tablename="knownGene",
                        miRBaseBuild="GRCh37")

## Mouse mm9
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="mm9",
                        tablename="knownGene")

## Mouse mm10
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="mm10",
                        tablename="ensGene")



## c elegans enseble trac
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="ce6",
                        tablename="ensGene")
                        
## Drosophila ensembl track
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="dm3",
                        tablename="ensGene")
                        
## Arabidopsis biomaRt
makeTxDbPackageFromBiomart(version=version,
                          maintainer=maintainer,
                          author=author,
		 	  destDir=TxDbOutDir,
                          biomart="plants_mart_12",
                          dataset="athaliana_eg_gene")
                                         

## Rat rn4 ensembl track
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="rn4",
                        tablename="ensGene")



## yeast ensembl genes is no longer available?
## Yeast sacCer2 ensembl track
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="sacCer2",
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
                                         


