### This function is called from makeTerminalDBPkgs.R which is part of the 
### release build code base. It creates new TxDb packages for a release and 
### updates existing packages that have modified tracks.

### Check for new or updated tracks at
###  https://genome.ucsc.edu/cgi-bin/hgGateway?hgsid=587155737_kHXGzURQSTDan47Ux2E8zNtqcI0A

TxDbPackagesForRelease <- 
    function(version, 
             destDir='.',
             maintainer= paste0("Bioconductor Package Maintainer ",
                         "<maintainer@bioconductor.org>"),
             author="Bioconductor Core Team")
{

    ##
    ## Start building packages for Bioconductor 3.7
    ##

    cat("building bosTau8\n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="bosTau8",
                            tablename="refGene")    

    cat("building ce11\n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="ce11",
                            tablename="refGene")    

    cat("building canFam3\n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="canFam3",
                            tablename="refGene")    

    cat("building danRer10\n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="danRer10",
                            tablename="refGene")    

    cat("building galGal4\n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="galGal4",
                            tablename="refGene")    

    cat("building galGal5\n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="galGal5",
                            tablename="refGene")    

    cat("building rheMac3\n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="rheMac3",
                            tablename="refGene")    

    cat("building rheMac8\n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="rheMac8",
                            tablename="refGene")    

    cat("building panTro4\n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="panTro4",
                            tablename="refGene")    

    cat("building rn5\n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="rn5",
                            tablename="refGene")    

    cat("building rn6\n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="rn6",
                            tablename="refGene")    

    cat("building susScr3\n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="susScr3",
                            tablename="refGene")    

    ##
    ## End building packages for Bioconductor 3.7
    ##
}
