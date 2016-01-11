###

if (FALSE) {
  ### Code use to obtain the hg19_knownGene_sample and
  ### hg19_knownToLocusLink_sample datasets.

  .get_hg19_knownGene <- function()
  {
    genome <- "hg19"
    tablename <- "knownGene"
    track <- GenomicFeatures:::.tablename2track(tablename, genome)
    library(rtracklayer)
    session <- browserSession()
    genome(session) <- genome
    track_tables <- tableNames(ucscTableQuery(session, track=track))
    query <- ucscTableQuery(session, track, table=tablename)
    getTable(query)
  }

  .get_hg19_knownToLocusLink <- function()
  {
    genome <- "hg19"
    tablename <- "knownGene"
    track <- GenomicFeatures:::.tablename2track(tablename, genome)
    library(rtracklayer)
    session <- browserSession()
    genome(session) <- genome
    mapdef <- GenomicFeatures:::.howToGetTxName2GeneIdMapping(tablename)
    txname2geneid <- GenomicFeatures:::.fetchTxName2GeneIdMappingFromUCSC(
                                 session, track, tablename, mapdef)
    txname2geneid$genes
  }

  library(GenomicFeatures)
  hg19_knownGene <- .get_hg19_knownGene()
  hg19_knownToLocusLink <- .get_hg19_knownToLocusLink()

  ## Pick up random samples.
  set.seed(333)
  sample_idx1 <- sort(sample(nrow(hg19_knownGene), 50))
  sample_idx2 <- which(hg19_knownToLocusLink$tx_name %in%
                       hg19_knownGene$name[sample_idx1])
  sample_idx2 <- which(hg19_knownToLocusLink$gene_id %in%
                       hg19_knownToLocusLink$gene_id[sample_idx2])
  hg19_knownToLocusLink_sample <-
                 droplevels(hg19_knownToLocusLink[sample_idx2, ])
  selected_tx_names <- union(hg19_knownGene$name[sample_idx1],
                             hg19_knownToLocusLink_sample$tx_name)
  sample_idx1 <- which(hg19_knownGene$name %in% selected_tx_names)
  hg19_knownGene_sample <- droplevels(hg19_knownGene[sample_idx1, ])

  ## Save them.
  save(hg19_knownToLocusLink_sample,
       file="GenomicFeatures/inst/extdata/hg19_knownToLocusLink_sample.rda")
  save(hg19_knownGene_sample,
       file="GenomicFeatures/inst/extdata/hg19_knownGene_sample.rda")
}

test_makeTxDbFromUCSCTxTable <- function()
{
    ## want
    txdb0_file <- system.file(
                      "extdata",
                      "hg19_knownGene_sample.sqlite",
                      package="GenomicFeatures"
                  )
    txdb0 <- loadDb(txdb0_file)

    ## get
    knownGene_sample_file <- system.file(
                                 "extdata",
                                 "hg19_knownGene_sample.rda",
                                 package="GenomicFeatures"
                             )
    load(knownGene_sample_file)
    knownToLocusLink_sample_file <- system.file(
                                        "extdata",
                                        "hg19_knownToLocusLink_sample.rda",
                                        package="GenomicFeatures"
                                    )
    load(knownToLocusLink_sample_file)
    genes <- data.frame(tx_name=hg19_knownToLocusLink_sample$tx_name,
                        gene_id=hg19_knownToLocusLink_sample$gene_id)
    txdb1 <- GenomicFeatures:::.makeTxDbFromUCSCTxTable(
                 hg19_knownGene_sample, genes,
                 "hg19", "knownGene", "UCSC Genes", "Entrez Gene ID", FALSE,
                 circ_seqs="chrM")

    ## compare
    ok <- GenomicFeatures:::compareTxDbs(txdb1, txdb0)
    checkTrue(ok)
}

