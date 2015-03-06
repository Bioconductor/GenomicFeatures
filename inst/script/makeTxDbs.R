library(GenomicFeatures)

## Best if I pass this down (from makeTerminal)
## version = "2.7.0" ## best if I pass this down (from makeTerminal)
## TxDbOutDir = paste(outDir,"_TxDbs",sep="")  

maintainer = "Bioconductor Package Maintainer <maintainer@bioconductor.org>"
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

## Human HG19 lincRNAs
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="hg19",
                        tablename="lincRNAsTranscripts",
                        miRBaseBuild="GRCh37")

## Human HG38
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="hg38",
                        tablename="knownGene")



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

## Mouse mm10 knownGene
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="mm10",
                        tablename="knownGene")


## c elegans enseble track
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
                          biomart="plants_mart_25",
                          dataset="athaliana_eg_gene")
                                         

## Rat rn4 ensembl track
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="rn4",
                        tablename="ensGene")

## Rat rn5 ensembl track
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="rn5",
                        tablename="refGene")



## yeast ensembl genes is no longer available?
## Yeast sacCer2 ensembl track
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="sacCer2",
                        tablename="sgdGene")

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
                                         


