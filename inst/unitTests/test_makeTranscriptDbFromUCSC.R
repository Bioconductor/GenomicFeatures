###

### The 2 utilities below were used to obtain the hg19_knownGene_sample
### and hg19_knownToLocusLink_sample datasets:
###   library(GenomicFeatures)
###   source(system.file("unitTests", "test_makeTranscriptDbFromUCSC.R",
###                      package="GenomicFeatures"))
###   hg19_knownGene_sample <- .make_hg19_knownGene_sample()
###   hg19_knownToLocusLink_sample <-
###           .make_hg19_knownToLocusLink_sample(hg19_knownGene_sample)
.make_hg19_knownGene_sample <- function()
{
    genome <- "hg19"
    tablename <- "knownGene"
    track <- GenomicFeatures:::.tablename2track(tablename, genome)
    library(rtracklayer)
    session <- browserSession()
    genome(session) <- genome
    track_tables <- tableNames(ucscTableQuery(session, track=track))
    query <- ucscTableQuery(session, track, table=tablename)
    hg19_knownGene <- getTable(query)
    set.seed(333)
    sample_idx <- sort(sample(nrow(hg19_knownGene), 100))
    droplevels(hg19_knownGene[sample_idx, ])
}

.make_hg19_knownToLocusLink_sample <- function(hg19_knownGene_sample)
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
    hg19_knownToLocusLink <- txname2geneid$genes
    sample_idx <- which(hg19_knownToLocusLink$tx_name %in%
                        hg19_knownGene_sample$name)
    droplevels(hg19_knownToLocusLink[sample_idx, ])
}

test_makeTranscriptDbFromUCSCTxTable <- function()
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
    txdb1 <- GenomicFeatures:::.makeTranscriptDbFromUCSCTxTable(
                 hg19_knownGene_sample, genes,
                 "hg19", "knownGene", "Entrez Gene ID", FALSE,
                 DEFAULT_CIRC_SEQS)

    ## compare
    ok <- GenomicFeatures:::compareTxDbs(txdb1, txdb0)
    checkTrue(ok)
}

