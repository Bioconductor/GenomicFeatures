test_transcripts <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))

    ## Test misuse.
    checkException(transcripts(data.frame()), silent=TRUE)
    checkException(transcripts(txdb, filter=list(bad=1:10)), silent=TRUE)
    checkException(transcripts(txdb, columns="bad"), silent=TRUE)

    ## Test 1.
    current <- transcripts(txdb, filter=list(gene_id="139231"))
    metadata(current) <- list()

    target <- GRanges("chrX",
                      ranges=IRanges(start=c(103411156, 103430747),
                                     end  =c(103440582, 103440582)),
                      strand=strand("+"),
                      tx_id=142:143,
                      tx_name=c("uc004elw.3", "uc004elx.3"),
                      seqinfo=seqinfo(txdb))
    checkIdentical(target, current)

    ## Test 2.
    filter <- list(tx_chrom=c("chr12", "chr14"), tx_strand="-")
    current <- transcripts(txdb, columns=c("tx_id", "tx_name",
                                           "exon_id", "exon_rank"),
                                 filter=filter)
    metadata(current) <- list()

    target <- GRanges("chr12",
                      ranges=IRanges(start=52753790, end=52761309),
                      strand=strand("-"),
                      tx_id=87L,
                      tx_name="uc001sag.3",
                      exon_id=IntegerList(334:326),
                      exon_rank=IntegerList(1:9),
                      seqinfo=seqinfo(txdb))
    checkIdentical(target, current)

    ## Test 3.
    filter <- list(gene_id=c("220004", "1183", "10186"))
    current <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"),
                                 filter=filter)
    metadata(current) <- list()

    target_tx_id <- c(91L, 136:137)
    target_gene_id <- CharacterList("10186", "1183", "1183")
    target <- GRanges(c("chr13", "chrX", "chrX"),
                      ranges=IRanges(start=c(39917029, 10124985, 10124985),
                                     end  =c(40177356, 10205699, 10205699)),
                      strand=strand(c("-", "+", "+")),
                      tx_id=target_tx_id,
                      tx_name=c("uc001uxf.3", "uc004csy.4", "uc011mid.3"),
                      gene_id=target_gene_id,
                      seqinfo=seqinfo(txdb))
    checkIdentical(target, current)
}

test_transcripts_after_seqlevelsStyle_switch <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                               package="GenomicFeatures"))
    checkIdentical(seqlevelsStyle(txdb), "UCSC")
    seqlevelsStyle(txdb) <- "NCBI"
    checkIdentical(seqlevelsStyle(txdb), c("NCBI", "UCSC"))

    filter <- list(gene_id=c("220004", "1183", "10186"))
    current <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"),
                                 filter=filter)
    metadata(current) <- list()

    checkIdentical(c("13", "X", "X"), as.character(seqnames(current)))

    target_tx_id <- c(91L, 136:137)
    target_gene_id <- CharacterList("10186", "1183", "1183")
    target <- GRanges(c("13", "X", "X"),
                      ranges=IRanges(start=c(39917029, 10124985, 10124985),
                                     end  =c(40177356, 10205699, 10205699)),
                      strand=strand(c("-", "+", "+")),
                      tx_id=target_tx_id,
                      tx_name=c("uc001uxf.3", "uc004csy.4", "uc011mid.3"),
                      gene_id=target_gene_id,
                      seqinfo=seqinfo(txdb))
    checkIdentical(target, current)    
}

test_exons <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                                      package="GenomicFeatures"))

    ## Test misuse.
    checkException(exons(data.frame()), silent=TRUE)
    checkException(exons(txdb, filter=list(bad=1:10)), silent=TRUE)

    ## Test 1.
    current <- exons(txdb, filter=list(tx_name="uc001gde.2"))
    metadata(current) <- list()

    target <- GRanges("chr1",
                      ranges=IRanges(start=c(165513478, 165532742),
                                     end  =c(165514155, 165533185)),
                      strand=strand("+"),
                      exon_id=29:30,
                      seqinfo=seqinfo(txdb))
    checkIdentical(target, current)

    ## Test 2.
    filter <- list(exon_chrom=c("chr5", "chr14"), exon_strand="-")
    current <- exons(txdb, columns=c("exon_id", "tx_name", "gene_id"),
                           filter=filter)
    metadata(current) <- list()

    target_ranges <-
        IRanges(start=c(134363424, 134366966, 134369403, 170732985),
                end  =c(134365011, 134367198, 134369964, 170735759))
    target <- GRanges("chr5",
                      ranges=target_ranges,
                      strand=strand("-"),
                      exon_id=182:185,
                      tx_name=CharacterList("uc010jea.3", "uc010jea.3",
                                            "uc010jea.3", "uc003mbe.2"),
                      gene_id=CharacterList("5307", "5307", "5307", NULL),
                      seqinfo=seqinfo(txdb))
    checkIdentical(target, current)
}

test_cds <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                         package="GenomicFeatures"))

    ## Test misuse.
    checkException(cds(data.frame()), silent=TRUE)
    checkException(cds(txdb, filter=list(bad=1:10)), silent=TRUE)

    ## Test 1.
    current <- cds(txdb, filter=list(tx_name="uc001gde.2"))
    metadata(current) <- list()

    target <- GRanges("chr1",
                      ranges=IRanges(start=c(165513534, 165532742),
                                     end  =c(165514155, 165533061)),
                      strand=strand("+"),
                      cds_id=23:24,
                      seqinfo=seqinfo(txdb))
    checkIdentical(target, current)

    ## Test 2.
    filter <- list(cds_chrom=c("chr5", "chr14"), cds_strand="-")
    current <- cds(txdb, columns=c("exon_id", "tx_name", "gene_id"),
                         filter=filter)
    metadata(current) <- list()

    target_ranges <-
        IRanges(start=c(134364469, 134366966, 134369403, 170735359),
                end  =c(134365011, 134367198, 134369571, 170735634))
    target <- GRanges("chr5",
                      ranges=target_ranges,
                      strand=strand("-"),
                      exon_id=as(182:185, "IntegerList"),
                      tx_name=CharacterList("uc010jea.3", "uc010jea.3",
                                            "uc010jea.3", "uc003mbe.2"),
                      gene_id=CharacterList("5307", "5307", "5307", NULL),
                      seqinfo=seqinfo(txdb))
    checkIdentical(target, current)
}

test_promoters <- function()
{
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                               package="GenomicFeatures"))
    tx <- transcripts(txdb, use.names=TRUE)
    current <- promoters(txdb)
    checkTrue(validObject(current))
    checkEquals(colnames(mcols(current)), c("tx_id", "tx_name"))
    checkIdentical(current, promoters(tx))
    current <- terminators(txdb)
    checkTrue(validObject(current))
    checkEquals(colnames(mcols(current)), c("tx_id", "tx_name"))
    checkIdentical(current, terminators(tx))
}


test_translateCols <- function(){
    txdb <- loadDb(system.file("extdata", "hg19_knownGene_sample.sqlite", 
                   package="GenomicFeatures"))
    tx1 <- transcripts(txdb, columns=c("tx_id", "tx_name", "cds_id"))
    checkEquals(colnames(mcols(tx1)), c("tx_id", "tx_name", "cds_id"))
    tx2 <- transcripts(txdb, columns=c("TXID", "TXNAME", "CDSID"))
    checkEquals(colnames(mcols(tx2)), c("TXID", "TXNAME", "CDSID"))
    tx3 <- transcripts(txdb, columns=c(bob="CDSID"))
    checkEquals(colnames(mcols(tx3)), c("bob"))
    tx4 <- transcripts(txdb, columns=c(bob="cds_id"))
    checkEquals(colnames(mcols(tx4)), c("bob"))
    ## And these two cases should both explode. ;)
    checkException(transcripts(txdb, columns=c("")))
    checkException(transcripts(txdb, columns=c("bob")))
}

