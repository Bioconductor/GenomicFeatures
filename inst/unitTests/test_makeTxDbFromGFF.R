## Test the helper functions, and then test the big one.

## 1st set up some resources for all the tests
## require(rtracklayer)

## require(GenomicFeatures)

## IF I require GenomicFeatures here it will 'fix' the problem for R -e
## 'BiocGenerics:::testPackage(pkgname="GenomicFeatures",
## pattern="test_makeTxDbFromGFF")' ...  but NOT for R CMD check
## GenomicFeatures_x.y.z.tar.gz
library(GenomicFeatures)
gffFile <- system.file("extdata","GFF3_files","a.gff3",package="GenomicFeatures")

gff3 <- rtracklayer:::import(gffFile, format="gff3", asRangedData=FALSE)
gff3 <- GenomicFeatures:::.makeSelectedFieldsLowerCase(gff3,
                               fieldsToCast=c('id','parent'))

gtfFile <- system.file("extdata","GTF_files","Aedes_aegypti.partial.gtf",
                       package="GenomicFeatures")

gtf <- rtracklayer:::import(gtfFile, format="gtf", asRangedData=FALSE)

flyFile <- system.file("extdata","GFF3_files","dmel-1000-r5.11.filtered.gff",
                       package="GenomicFeatures")

## bad bacterial GFFs require use of special argument to ignore most of data.
gffB <- system.file("extdata","GFF3_files","NC_011025.gff",package="GenomicFeatures")


## test .deduceExonRankings (usually used by gff) 
test_deduceExonRankings <- function(){
  data <- GenomicFeatures:::.prepareGFF3data.frame(gff3,
                                                   exonRankAttributeName=NA,
                                                   gffGeneIdAttributeName=NA)
  exs <- GenomicFeatures:::.prepareGFF3Fragments(data,type="exon")
  suppressWarnings(
    res <- GenomicFeatures:::.deduceExonRankings(exs))
  ## TESTING
  checkTrue(!any(is.na(res$exon_rank))) ## is filled out.
  checkTrue(is.integer(res$exon_rank)) ## is integers
  checkIdentical(res$exon_rank[3:5],1:3)  ## correct val for + strand set
  checkIdentical(res$exon_rank[10:14],5:1)  ## correct val for - strand set
}


## test .mergeFramesViaRanges (used by both)
test_mergeFramesViaRanges <- function(){
  data <- GenomicFeatures:::.prepareGFF3data.frame(gff3,
                                                   exonRankAttributeName=NA,
                                                   gffGeneIdAttributeName=NA)

  exs <- GenomicFeatures:::.prepareGFF3Fragments(data,type="exon")
  cds <- GenomicFeatures:::.prepareGFF3Fragments(data,type="CDS")
    suppressWarnings(
      exs <- GenomicFeatures:::.deduceExonRankings(exs))
  res <- GenomicFeatures:::.mergeFramesViaRanges(exs, cds)
  ## TESTING
  checkTrue(!any(is.na(res$exon_rank))) ## rank info is filled out.
  expectedCols <- c("exon_chrom","exon_start","exon_end","exon_strand","type",
    "exon_name","tx_name","exon_rank","cds_chrom","cds_start",
    "cds_end","cds_strand","type","cds_name","tx_name","exon_rank")
  checkIdentical(expectedCols, colnames(res)) ## expected fields are present
  ## after matching, cds should always be 'within' the exons...
  checkTrue(all(res$exon_start <= res$cds_start, na.rm=TRUE)) 
  checkTrue(all(res$exon_end >= res$cds_end, na.rm=TRUE))

  ## now test our gtf file
  data <- GenomicFeatures:::.prepareGTFdata.frame(gtf,
                                           exonRankAttributeName="exon_number") 
  exs <- GenomicFeatures:::.prepareGTFFragments(data,type="exon")
  cds <- GenomicFeatures:::.prepareGTFFragments(data,type="CDS")
  res <- GenomicFeatures:::.mergeFramesViaRanges(exs, cds)
  ## TESTING
  checkTrue(!any(is.na(res$exon_rank))) ## rank info is filled out.
  expectedCols <- c("tx_name","exon_rank","exon_chrom","exon_strand",
    "exon_start","exon_end","tx_name","exon_rank","cds_chrom","cds_strand",
    "cds_start","cds_end")
  checkIdentical(expectedCols, colnames(res)) ## expected fields are present
  ## after matching, cds should always be 'within' the exons...
  checkTrue(all(res$exon_start <= res$cds_start, na.rm=TRUE))
  checkTrue(all(res$exon_end >= res$cds_end, na.rm=TRUE))
}


## test .deduceTranscriptsFromGTF (used by gtf)
test_deduceTranscriptsFromGTF <- function(){
  data <- GenomicFeatures:::.prepareGTFdata.frame(gtf,
                                           exonRankAttributeName="exon_number")
  suppressWarnings(
   res <- GenomicFeatures:::.deduceTranscriptsFromGTF(data))
  ## TESTING
  checkTrue(!any(is.na(res$tx_start))) ## rank info is filled out.
  checkTrue(!any(is.na(res$tx_end))) ## rank info is filled out.
  expectedCols <- c("tx_id","tx_name","tx_chrom","tx_strand","tx_start",
                    "tx_end","gene_id")
  checkIdentical(expectedCols, colnames(res)) ## expected fields are present
  checkTrue(all(res$tx_end > res$tx_start, na.rm=TRUE))
}



## test  .prepareGTFTables
test_prepareGTFTables <- function(){
  suppressWarnings(
  res <- GenomicFeatures:::.prepareGTFTables(gtf,
                                          exonRankAttributeName="exon_number"))
  ## TESTING
  checkTrue(length(res)==3)
  checkTrue(class(res$transcripts)=="data.frame")
  checkTrue(dim(res$transcripts)[1] == 105)
  checkTrue(dim(res$transcripts)[2] == 6)
  checkTrue(class(res$genes)=="data.frame")
  checkTrue(dim(res$genes)[1] == 105)
  checkTrue(dim(res$genes)[2] == 2)
  checkTrue(class(res$splicings)=="data.frame")
  checkTrue(dim(res$splicings)[1] == 414)
  checkTrue(dim(res$splicings)[2] == 8)
}



## test  .prepareGFF3Tables
test_prepareGFF3Tables <- function(){
  suppressWarnings(
  res <- GenomicFeatures:::.prepareGFF3Tables(gff3,
                                              exonRankAttributeName=NA,
                                              gffGeneIdAttributeName=NA,
                                              useGenesAsTranscripts=FALSE))
  ## TESTING
  checkTrue(length(res)==3)
  checkTrue(class(res$transcripts)=="data.frame")
  checkTrue(dim(res$transcripts)[1] == 488)
  checkTrue(dim(res$transcripts)[2] == 6)
  checkTrue(class(res$genes)=="data.frame")
  checkTrue(dim(res$genes)[1] == 488)
  checkTrue(dim(res$genes)[2] == 2)
  checkTrue(class(res$splicings)=="data.frame")
  checkTrue(dim(res$splicings)[1] == 1268)
  ## only fails when run as a unit test.
  checkTrue(dim(res$splicings)[2] == 9)
}




## Test that outputs match what is expected. ## BOOM
test_makeTxDbFromGFF <- function(){  
  ## wanted
  gffDBFile <- system.file("extdata","TESTGFF.sqlite",
                           package="GenomicFeatures")
  txdb_gff <- loadDb(gffDBFile)

  ## generated
  suppressWarnings(
  txdb <- makeTxDbFromGFF(file=gffFile,
         format="gff3",
         dataSource="partial GFF file for Tomatoes for testing",
         species="Solanum lycopersicum")
   )

  ## test
  checkTrue(GenomicFeatures:::compareTxDbs(txdb, txdb_gff))

  
  ## wanted
  gtfDBFile <- system.file("extdata","TESTGTF.sqlite",
                           package="GenomicFeatures")
  txdb_gtf <- loadDb(gtfDBFile)

  ## generated
  chrominfo <- data.frame(chrom = c('supercont1.1','supercont1.2'),
                        length=c(5220442, 5300000),
                        is_circular=c(FALSE, FALSE))
  
  suppressWarnings(
  txdb2 <- makeTxDbFromGFF(file=gtfFile,
         format="gtf",
         exonRankAttributeName="exon_number",
         chrominfo= chrominfo,
         dataSource=paste("ftp://ftp.ensemblgenomes.org/pub/metazoa/",
                          "release-13/gtf/aedes_aegypti/",sep=""),
         species="Aedes aegypti"))

  ## test
  checkTrue(GenomicFeatures:::compareTxDbs(txdb2, txdb_gtf))


  ## wanted
  flyDBFile <- system.file("extdata","TESTFLYGFF3.sqlite",
                           package="GenomicFeatures")
  txdb_fly <- loadDb(flyDBFile)

  suppressWarnings(
  txdb3 <- makeTxDbFromGFF(file=flyFile,
                                   format="gff3",
                                   dataSource="gff file from flybase",
                                   gffGeneIdAttributeName = "geneID",
                                   species="Drosophila melanogaster")
                   )
  
  ## test
  ## print(txdb3)
  ## print(lst1 <- lapply(as.list(txdb3), head, n=30))
  ## print(lst2 <- lapply(as.list(txdb_fly), head, n=30))
  ## print(all.equal(lst1, lst2))
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
                                      species = "Mycoplasma arthritidis",
                                      useGenesAsTranscripts=TRUE)

  ## Tests
  checkTrue(class(txdb_bac)=="TxDb")
  checkTrue(length(transcripts(txdb_bac)) > 100)

}



## the gffGeneIdAttributeName argument does not seem to work?  Or how does it work?  TODO: add some tests for this.  Also: it looks like it might be a little different from what I was about to add (which would be values put into the genes row of the genes records (Name or dbxRef)  (this just means that I probably should allow other argument to indicate which type of record to search in).

