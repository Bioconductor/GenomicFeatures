### This function is called from makeTerminalDBPkgs.R which is part of the 
### release build code base. It creates new TxDb packages for a release and 
### updates existing packages that have modified tracks.

TxDbPackagesForRelease <- 
    function(version, 
             destDir='.',
             maintainer= paste0("Bioconductor Package Maintainer ",
                         "<maintainer@bioconductor.org>"),
             author="Bioconductor Core Team")
{
    ## New packages for Bioconductor 3.3
    cat("building galGal4 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=txdbDir,
                            genome="galGal4",
                            tablename="refGene")
    cat("building ce11 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=txdbDir,
                            genome="ce11",
                            tablename="refGene")
    cat("building dm6 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=txdbDir,
                            genome="dm6",
                            tablename="ensGene")

    cat("building rn6 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=txdbDir,
                            genome="rn6",
                            tablename="refGene")

    cat("building panTro4 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=txdbDir,
                            genome="panTro4",
                            tablename="refGene")

    cat("building bostau8 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=txdbDir,
                            genome="bostau8",
                            tablename="refGene")

    cat("building canFam3 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=txdbDir,
                            genome="canFam3",
                            tablename="refGene")

    cat("building danRer10 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=txdbDir,
                            genome="danRer10",
                            tablename="refGene")
 
    cat("building rheMac3 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=txdbDir,
                            genome="rheMac3",
                            tablename="refGene")

    cat("building susScr3 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=txdbDir,
                            genome="susScr3",
                            tablename="refGene")

    ## Update existing packages
    cat("building rn5 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=txdbDir,
                            genome="rn5",
                            tablename="refGene")
}

### Static packages - no updates needed.
#
### Human HG18
#cat("building hg18 \n")
#makeTxDbPackageFromUCSC(version=version,
#                        maintainer=maintainer,
#                        author=author,
#			destDir=txdbDir,
#                        genome="hg18",
#                        tablename="knownGene")
#
### Human HG19
#cat("building hg19 \n")
#makeTxDbPackageFromUCSC(version=version,
#                        maintainer=maintainer,
#                        author=author,
#			destDir=txdbDir,
#                        genome="hg19",
#                        tablename="knownGene",
#                        miRBaseBuild="GRCh37")
#
### Human HG19 lincRNAs
#cat("building hg19 lincRNAs \n")
#makeTxDbPackageFromUCSC(version=version,
#                        maintainer=maintainer,
#                        author=author,
#			destDir=txdbDir,
#                        genome="hg19",
#                        tablename="lincRNAsTranscripts",
#                        miRBaseBuild="GRCh37")
#
### Human HG38
#makeTxDbPackageFromUCSC(version=version,
#                        maintainer=maintainer,
#                        author=author,
#			destDir=txdbDir,
#                        genome="hg38",
#                        tablename="knownGene")
#
### Mouse mm9
#cat("building mm9 \n")
#makeTxDbPackageFromUCSC(version=version,
#                        maintainer=maintainer,
#                        author=author,
#			destDir=txdbDir,
#                        genome="mm9",
#                        tablename="knownGene")
#
### Mouse mm10
#cat("building mm10 \n")
#makeTxDbPackageFromUCSC(version=version,
#                        maintainer=maintainer,
#                        author=author,
#			destDir=txdbDir,
#                        genome="mm10",
#                        tablename="ensGene")
#
### Mouse mm10 knownGene
#cat("building mm10 knownGene \n")
#makeTxDbPackageFromUCSC(version=version,
#                        maintainer=maintainer,
#                        author=author,
#			destDir=txdbDir,
#                        genome="mm10",
#                        tablename="knownGene")
#
### c elegans enseble track
#cat("building ce6 \n")
#makeTxDbPackageFromUCSC(version=version,
#                        maintainer=maintainer,
#                        author=author,
#			destDir=txdbDir,
#                        genome="ce6",
#                        tablename="ensGene")
# 
### Drosophila ensembl track
#cat("building dm3 \n")
#makeTxDbPackageFromUCSC(version=version,
#                        maintainer=maintainer,
#                        author=author,
#			destDir=txdbDir,
#                        genome="dm3",
#                        tablename="ensGene")
#
### Arabidopsis biomaRt
#cat("building plants_mart \n")
#makeTxDbPackageFromBiomart(version=version,
#                          maintainer=maintainer,
#                          author=author,
#		 	  destDir=txdbDir,
#                          host="plants.ensembl.org",
#                          biomart="plants_mart",
#                          dataset="athaliana_eg_gene")
#
### Rat rn4 ensembl track
#cat("building rn4 \n")
#makeTxDbPackageFromUCSC(version=version,
#                        maintainer=maintainer,
#                        author=author,
#			destDir=txdbDir,
#                        genome="rn4",
#                        tablename="ensGene")
#
### Yeast sacCer2 ensembl track
#cat("building sacCer2 sgdGene \n")
#makeTxDbPackageFromUCSC(version=version,
#                        maintainer=maintainer,
#                        author=author,
#			destDir=txdbDir,
#                        genome="sacCer2",
#                        tablename="sgdGene")
#
#cat("building sacCer3 \n")
#makeTxDbPackageFromUCSC(version=version,
#                        maintainer=maintainer,
#                        author=author,
#			destDir=txdbDir,
#                        genome="sacCer3",
#                        tablename="sgdGene")
#
#cat("building sacCer2 ensGene \n")
#makeTxDbPackageFromUCSC(version=version,
#                        maintainer=maintainer,
#                        author=author,
#                        destDir=txdbDir,
#                        genome="sacCer2",
#                        tablename="ensGene")
