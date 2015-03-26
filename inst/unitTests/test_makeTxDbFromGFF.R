## IF I require GenomicFeatures here it will 'fix' the problem for R -e
## 'BiocGenerics:::testPackage(pkgname="GenomicFeatures",
## pattern="test_makeTxDbFromGFF")' ...  but NOT for R CMD check
## GenomicFeatures_x.y.z.tar.gz
library(GenomicFeatures)
gffFile <- system.file("extdata","GFF3_files","a.gff3",package="GenomicFeatures")

gtfFile <- system.file("extdata","GTF_files","Aedes_aegypti.partial.gtf",
                       package="GenomicFeatures")

flyFile <- system.file("extdata","GFF3_files","dmel-1000-r5.11.filtered.gff",
                       package="GenomicFeatures")

## bad bacterial GFFs require use of special argument to ignore most of data.
gffB <- system.file("extdata","GFF3_files","NC_011025.gff",package="GenomicFeatures")

## Test that outputs match what is expected. ## BOOM
test_makeTxDbFromGFF <- function(){  
  ## wanted
  gffDBFile <- system.file("extdata", "GFF3_files", "a.sqlite",
                           package="GenomicFeatures")
  txdb_gff <- loadDb(gffDBFile)

  ## generated
  txdb <- makeTxDbFromGFF(file=gffFile,
         format="gff3",
         dataSource="partial GFF file for Tomatoes for testing",
         organism="Solanum lycopersicum")

  ## test
  checkTrue(GenomicFeatures:::compareTxDbs(txdb, txdb_gff))

  
  ## wanted
  gtfDBFile <- system.file("extdata", "GTF_files",
                           "Aedes_aegypti.partial.sqlite",
                           package="GenomicFeatures")
  txdb_gtf <- loadDb(gtfDBFile)

  ## generated
  chrominfo <- data.frame(chrom = c('supercont1.1','supercont1.2'),
                        length=c(5220442, 5300000),
                        is_circular=c(FALSE, FALSE))
  
  txdb2 <- makeTxDbFromGFF(file=gtfFile,
         format="gtf",
         chrominfo= chrominfo,
         dataSource=paste("ftp://ftp.ensemblgenomes.org/pub/metazoa/",
                          "release-13/gtf/aedes_aegypti/",sep=""),
         organism="Aedes aegypti")

  ## test
  checkTrue(GenomicFeatures:::compareTxDbs(txdb2, txdb_gtf))


  ## wanted
  flyDBFile <- system.file("extdata", "GFF3_files",
                           "dmel-1000-r5.11.filtered.sqlite",
                           package="GenomicFeatures")
  txdb_fly <- loadDb(flyDBFile)

  txdb3 <- makeTxDbFromGFF(file=flyFile,
                                   format="gff3",
                                   dataSource="gff file from flybase",
                                   organism="Drosophila melanogaster")
  
  checkTrue(GenomicFeatures:::compareTxDbs(txdb3, txdb_fly))


  ## test for broken NCBI bacterial GFFs (that only seem to have
  ## reliable gene info and little else)
  chrominfoBac <- data.frame(chrom = c('NC_011025.1'),
                          length=c(830000), ## placeholder = iow it big enough
                          is_circular=c(TRUE))

  ## mostly I want to see if if can run this:
  txdb_bac <- makeTxDbFromGFF(file = gffB, format = "gff",
                                      chrominfo = chrominfoBac,
                                      dataSource = "NCBI",
                                      organism = "Mycoplasma arthritidis")

  ## Tests
  checkTrue(class(txdb_bac)=="TxDb")
  checkTrue(length(transcripts(txdb_bac)) > 100)
}

