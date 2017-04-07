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
    ## Build new tracks for Bioconductor 3.5
    cat("building galGal5 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="galGal5",
                            tablename="refGene")


    ## Update live tracks for Bioconductor 3.5
    cat("building ce11 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="ce11",
                            tablename="ensGene")
    cat("building dm6 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="dm6",
                            tablename="ensGene")
    cat("building rn5 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="rn5",
                            tablename="refGene")

    cat("building rn6 \n")
    makeTxDbPackageFromUCSC(version=version,
                            maintainer=maintainer,
                            author=author,
                            destDir=destDir,
                            genome="rn6",
                            tablename="refGene")
}
