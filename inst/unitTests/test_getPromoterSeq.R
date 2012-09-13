# GenomicFeatures/inst/unitTests/test_getPromoterSeqs.R
#--------------------------------------------------------------------------------
library (RUnit)
library (BSgenome.Hsapiens.UCSC.hg19)
library (TxDb.Hsapiens.UCSC.hg19.knownGene)
#--------------------------------------------------------------------------------
e2f3 <- '1871'   # a gene on the plus strand
grb2 <- '2885'   # a gene on the minus strand
#--------------------------------------------------------------------------------
run.tests = function ()
{
  test.GRangesList.BSgenome.getPromoterSeq ()
  test.GRanges.BSgenome.getPromoterSeq ()
  test.GRangesList.Fasta.getPromoterSeq ()
  test.GRanges.Fasta.getPromoterSeq ()

} # run.tests
#--------------------------------------------------------------------------------
test.GRangesList.BSgenome.getPromoterSeq = function ()
{
  print ('--- test.GRangesList.BSgenome.getPromoterSeq')

  genes <- c (e2f3, grb2)
  transcriptCoordsByGene.GRangesList <- transcriptsBy (TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene") [genes]
  checkEquals (names (transcriptCoordsByGene.GRangesList), genes)

  promoter.seqs <- getPromoterSeq (transcriptCoordsByGene.GRangesList, Hsapiens, upstream=10, downstream=0)
  checkTrue (is (promoter.seqs, 'DNAStringSetList'))
  checkEquals (length (promoter.seqs), 2)
  checkEquals (names (promoter.seqs), genes)

  checkEquals (width (unlist(promoter.seqs)), rep (10, 5))
  checkEquals (as.character(unlist(promoter.seqs, use.names=F)),
              c ("GCTTCCTGGA", "GCTTCCTGGA", "CGGAGCCAGG", "CCTCGTGGAG", "CCTCGTGGAG"))

} # test.GRangesList.BSgenome.getPromoterSeq 
#--------------------------------------------------------------------------------
test.GRangesList.Fasta.getPromoterSeq = function ()
{
  print ('--- test.GRangesList.Fasta.getPromoterSeq')

  filename = file.path (find.package ('GenomicFeatures'), 'extdata', 'chr17.fa')
  checkTrue (file.exists (filename))
  sequence.from.fasta = open (Rsamtools::FaFile (filename))

  armc7 <- '79637'   # a gene on the plus strand, chr17:70,617,676-70,637,955
  grb2  <- '2885'   # a gene on the minus strand,  chr17:70,825,751 70,913,385

  genes <- c (armc7, grb2)
  transcriptCoordsByGene.GRangesList <- transcriptsBy (TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene") [genes]
  checkEquals (names (transcriptCoordsByGene.GRangesList), genes)

  promoter.seqs <- getPromoterSeq (transcriptCoordsByGene.GRangesList, sequence.from.fasta, upstream=10, downstream=10)
  checkTrue (is (promoter.seqs, 'DNAStringSetList'))
  checkEquals (length (promoter.seqs), 2)
  checkEquals (names (promoter.seqs), genes)
  
  checkEquals (width (unlist(promoter.seqs)), rep (20, 5))
  
  checkEquals (as.character (as.character(unlist(promoter.seqs, use.names=F))),
               c ('GCCCAGGAGCCGGGGACGGT', 'GCCCAGGAGCCGGGGACGGT', 'CGGAAGTCACAGAGCAGTCC',
                  'CCTCGTGGAGAAGTTCTCGC', 'CCTCGTGGAGAAGTTCTCGC'))

} # test.GRangesList.Fasta.getPromoterSeq 
#--------------------------------------------------------------------------------
test.GRanges.Fasta.getPromoterSeq = function ()
{
  print ('--- test.GRanges.Fasta.getPromoterSeq')

  filename = file.path (find.package ('GenomicFeatures'), 'extdata', 'chr17.fa')
  checkTrue (file.exists (filename))
  sequence.from.fasta = open (Rsamtools::FaFile (filename))

  transcriptCoordsByGene.GRanges <- transcriptsBy (TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene") [[grb2]]
  checkTrue(is (transcriptCoordsByGene.GRanges, 'GRanges'))
  checkTrue (is.null (names (transcriptCoordsByGene.GRanges)))  # wd have names only if its a list
  checkEquals (dim (mcols (transcriptCoordsByGene.GRanges)), c (2, 2))
  checkEquals (colnames (mcols (transcriptCoordsByGene.GRanges)), c ('tx_id', 'tx_name'))

  promoter.seqs <- getPromoterSeq (transcriptCoordsByGene.GRanges, sequence.from.fasta, upstream=10, downstream=10)
  checkTrue (is (promoter.seqs, 'DNAStringSet'))
  checkEquals (length (promoter.seqs), 2)

   # todo: the BSgenome version has names which are NA.  here we get the chromosome names.  why the difference?
   # checkTrue (is.null (names (promoter.seqs)))

  checkEquals (width (promoter.seqs), rep (20, 2))
  checkEquals (as.character (as.character (promoter.seqs)), c ('CCTCGTGGAGAAGTTCTCGC', 'CCTCGTGGAGAAGTTCTCGC'))
    # should be one more column in the metadata than in the metadata 
  checkEquals (dim (mcols (promoter.seqs)), c (2, 3))
  checkEquals (colnames (mcols (promoter.seqs)), c ('tx_id', 'tx_name', 'geneID'))

     # the input, a GRanges, had no names -- which are the source of geneID when the GRangesList version of
     # this methods is called.  so ensure that this lack of information was passed along into the
     # metadata of the returned promoter.seqs
  checkTrue (all (is.na (mcols(promoter.seqs)$geneID)))

} # test.GRanges.Fasta.getPromoterSeq 
#--------------------------------------------------------------------------------
test.GRanges.BSgenome.getPromoterSeq = function ()
{
  print ('--- test.GRanges.BSgenome.getPromoterSeq')

  transcriptCoordsByGene.GRanges <- transcriptsBy (TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene") [[e2f3]]
  checkTrue(is (transcriptCoordsByGene.GRanges, 'GRanges'))
  checkTrue (is.null (names (transcriptCoordsByGene.GRanges)))  # wd have names only if its a list
  checkEquals (dim (mcols (transcriptCoordsByGene.GRanges)), c (3, 2))
  checkEquals (colnames (mcols (transcriptCoordsByGene.GRanges)), c ('tx_id', 'tx_name'))

  promoter.seqs <- getPromoterSeq (transcriptCoordsByGene.GRanges, Hsapiens, upstream=10, downstream=0)
  checkTrue (is (promoter.seqs, 'DNAStringSet'))
  checkEquals (length (promoter.seqs), 3)
  checkTrue (is.null (names (promoter.seqs)))
  checkEquals (width (promoter.seqs), rep (10, 3))
  checkEquals (as.character (promoter.seqs), c ("GCTTCCTGGA", "GCTTCCTGGA", "CGGAGCCAGG"))

    # should be one more column in the metadata than in the metadata 
  checkEquals (dim (mcols (promoter.seqs)), c (3, 3))
  checkEquals (colnames (mcols (promoter.seqs)), c ('tx_id', 'tx_name', 'geneID'))

     # the input, a GRanges, had no names -- which are the source of geneID when the GRangesList version of
     # this methods is called.  so ensure that this lack of information was passed along into the
     # metadata of the returned promoter.seqs
  checkTrue (all (is.na (mcols(promoter.seqs)$geneID)))

} # test.GRanges.BSgenome.getPromoterSeq 
#--------------------------------------------------------------------------------
