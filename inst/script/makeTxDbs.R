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
cat("building plants_mart \n")
makeTxDbPackageFromBiomart(version=version,
                          maintainer=maintainer,
                          author=author,
		 	  destDir=TxDbOutDir,
                          host="plants.ensembl.org",
                          biomart="plants_mart",
                          dataset="athaliana_eg_gene")
                        
#Warning message:
#In .warning_on_BioMart_data_anomaly(bm_result, idx, id_prefix, msg) :
#  BioMart data anomaly: in the following transcripts, 
#  the CDS/UTR genomic coordinates are inconsistent with the 
#  "cds_start" and "cds_end" attributes.
#  (Showing only the first 3 out of 40 transcripts.)
                 

## Rat rn4 ensembl track
cat("building rn4 \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="rn4",
                        tablename="ensGene")

#Warning message:
#In .extractCdsLocsFromUCSCTxTable(ucsc_txtable, exon_locs) :
#  UCSC data anomaly in 71 transcript(s): the cds cumulative length is not
#  a multiple of 3 for transcripts ‘ENSRNOT00000064456’
#  ‘ENSRNOT00000056007’ ‘ENSRNOT00000055902’ ‘ENSRNOT00000001880’


## Rat rn5 ensembl track
cat("building rn5 \n")
makeTxDbPackageFromUCSC(version=version,
                        maintainer=maintainer,
                        author=author,
			destDir=TxDbOutDir,
                        genome="rn5",
                        tablename="refGene")


#Warning message:
#In .extractCdsLocsFromUCSCTxTable(ucsc_txtable, exon_locs) :
#  UCSC data anomaly in 1179 transcript(s): the cds cumulative length is
#  not a multiple of 3 for transcripts ‘NM_031007’ ‘NM_001113372’
#  ‘NM_001077677’ ‘NM_001113371’ ‘NM_032071’ ‘NM_012693’ ‘NM_031069’
#


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
