### =========================================================================
### makeTranscriptDbFromBiomart()
### -------------------------------------------------------------------------
###
### For people who want to tap BioMart.
### Typical use:
###   txdb <- makeTranscriptDbFromBiomart(biomart="ensembl",
###                                       dataset="hsapiens_gene_ensembl")
### Speed:
###   - for biomart="ensembl" and dataset="hsapiens_gene_ensembl":
###       (1) download takes about 8 min.
###       (2) db creation takes about 60-65 sec.
###

.extractCdsRangesFromBiomartTable <- function(bm_table)
{
    strand <- bm_table[["strand"]]
    cds_start <- exon_start <- bm_table[["exon_chrom_start"]]
    cds_end <- exon_end <- bm_table[["exon_chrom_end"]]
    utr5_start <- bm_table[["5_utr_start"]]
    utr5_end <- bm_table[["5_utr_end"]]
    utr3_start <- bm_table[["3_utr_start"]]
    utr3_end <- bm_table[["3_utr_end"]]

    if (!all(strand %in% c(1, -1)))
        stop("BioMart data anomaly: \"strand\" attribute should be 1 or -1")
    if (!is.numeric(exon_start)
     || !is.numeric(exon_end)
     || !is.numeric(utr5_start)
     || !is.numeric(utr5_end)
     || !is.numeric(utr3_start)
     || !is.numeric(utr3_end))
        stop("BioMart data anomaly: exon or utr coordinates don't ",
             "have a numeric type")
    no_utr5 <- is.na(utr5_start)
    if (!identical(no_utr5, is.na(utr5_end)))
        stop("BioMart data anomaly: NAs in \"5_utr_start\" attribute ",
             "don't match NAs in \"5_utr_end\" attribute")
    if (!all(utr5_start <= utr5_end, na.rm=TRUE))
        stop("BioMart data anomaly: some 5' UTR have a start > end")
    if (!all(utr5_start >= exon_start, na.rm=TRUE)
     || !all(utr5_end <= exon_end, na.rm=TRUE))
        stop("BioMart data anomaly: some 5' UTR are not within the exon limits")
    no_utr3 <- is.na(utr3_start)
    if (!identical(no_utr3, is.na(utr3_end)))
        stop("BioMart data anomaly: NAs in \"3_utr_start\" attribute ",
             "don't match NAs in \"3_utr_end\" attribute")
    if (!all(utr3_start <= utr3_end, na.rm=TRUE))
        stop("BioMart data anomaly: some 3' UTR have a start > end")
    if (!all(utr3_start >= exon_start, na.rm=TRUE)
     || !all(utr3_end <= exon_end, na.rm=TRUE))
        stop("BioMart data anomaly: some 3' UTR are not within the exon limits")

    idx <- strand == 1 & !no_utr5
    if (!all(utr5_start[idx] == exon_start[idx]))
        stop("BioMart data anomaly: some 5' UTR on the plus strand ",
             "don't start where the exon starts")
    cds_start[idx] <- utr5_end[idx] + 1L
    idx <- strand == 1 & !no_utr3
    if (!all(utr3_end[idx] == exon_end[idx]))
        stop("BioMart data anomaly: some 3' UTR on the plus strand ",
             "don't end where the exon ends")
    cds_end[idx] <- utr3_start[idx] - 1L
    idx <- strand == -1 & !no_utr3
    if (!all(utr3_start[idx] == exon_start[idx]))
        stop("BioMart data anomaly: some 3' UTR on the minus strand ",
             "don't start where the exon starts")
    cds_start[idx] <- utr3_end[idx] + 1L
    idx <- strand == -1 & !no_utr5
    if (!all(utr5_end[idx] == exon_end[idx]))
        stop("BioMart data anomaly: some 5' UTR on the minus strand ",
             "don't end where the exon ends")
    cds_end[idx] <- utr5_start[idx] - 1L
    ans <- IRanges(start=cds_start, end=cds_end)
    if (length(ans) != 0L) {
        cds_cumlength <-
            sapply(split(width(ans), bm_table$ensembl_transcript_id), sum)
        if (!all(cds_cumlength[as.vector(bm_table$ensembl_transcript_id)]
                 == bm_table$cds_length, na.rm=TRUE))
            stop("BioMart data anomaly: for some transcripts, the cds ",
                 "cumulative length inferred from the exon and UTR info ",
                 "doesn't match the \"cds_length\" attribute from BioMart")
        #if (!all(cds_cumlength %% 3L == 0L))
        #    warning("BioMart data anomaly: for some transcripts, the cds ",
        #            "cumulative length (\"cds_length\" attribute) is not ",
        #            "a multiple of 3")
    }
    ans
}

### Here is the status of some of the available BioMart databases returned by
### listMarts() as of 2/18/2010:
###
### biomart: ensembl
### version: ENSEMBL 56 GENES (SANGER UK)
###   50 datasets, one per organism, all named *_gene_ensembl.
###   Examples: hsapiens_gene_ensembl, mmusculus_gene_ensembl,
###             btaurus_gene_ensembl, cfamiliaris_gene_ensembl,
###             drerio_gene_ensembl, ggallus_gene_ensembl,
###             dmelanogaster_gene_ensembl, celegans_gene_ensembl,
###             scerevisiae_gene_ensembl, ...
###   -> they are ALL GOOD
###
### biomart: bacterial_mart_3
### version: ENSEMBL BACTERIA 3 (EBI UK)
###   134 datasets, all named *_gene
###   -> they all lack CDS/UTR information
###
### biomart: fungal_mart_3
### version: ENSEMBL FUNGAL 3 (EBI UK)
###   10 datasets, all named *_eg_gene
###   -> they all lack CDS/UTR information
###
### biomart: metazoa_mart_3
### version: ENSEMBL METAZOA 3 (EBI UK)
###   21 datasets, all named *_eg_gene
###   -> they all lack CDS/UTR information
###
### biomart: plant_mart_3
### version: ENSEMBL PLANT 3 (EBI UK)
###   Arabidopsis is here ("athaliana_eg_gene" dataset)
###   8 datasets, all named *_eg_gene
###   -> they all lack CDS/UTR information
###
### biomart: protist_mart_3
### version: ENSEMBL PROTISTS 3 (EBI UK)
###   3 datasets, all named *_gene
###   -> they all lack CDS/UTR information
###
### biomart: wormbase_current
### version: WORMBASE (CSHL US)
###   ?? (TODO)
###
### biomart: ENSEMBL_MART_ENSEMBL
### version: GRAMENE 30 ENSEMBL GENES
###   ?? (TODO)
###
### Note that listMarts() and listDatasets() are returning data frames where
### the columns are character factors for the former and "AsIs" character
### vectors for the latter.
makeTranscriptDbFromBiomart <- function(biomart="ensembl",
                                        dataset="hsapiens_gene_ensembl",
                                        transcript_ids=NULL)
{
    ## Could be that the user got the 'biomart' and/or 'dataset' values
    ## programmatically via calls to listMarts() and/or listDatasets().
    if (is.factor(biomart))
        biomart <- as.character(biomart)
    if (is(dataset, "AsIs"))
        dataset <- as.character(dataset)
    if (!isSingleString(biomart))
        stop("'biomart' must be a single string")
    if (biomart != "ensembl")
        stop("only 'biomart=\"ensembl\"' is supported for now")
    mart <- useMart(biomart=biomart, dataset=dataset)
    datasets <- listDatasets(mart)
    dataset_row <- which(as.character(datasets$dataset) == dataset)
    ## This should never happen (the above call to useMart() would have failed
    ## in the first place).
    if (length(dataset_row) != 1L)
        stop("the BioMart database \"", biomart, "\" has 0 (or ",
             "more than 1) \"", dataset, "\" datasets")
    description <- as.character(datasets$description)[dataset_row]
    version <- as.character(datasets$version)[dataset_row]
    if (is.null(transcript_ids)) {
        filters <- values <- ""
    } else if (is.character(transcript_ids)
            && !any(is.na(transcript_ids))) {
        filters <- "ensembl_transcript_id"
        values <- transcript_ids
    } else {
            stop("'transcript_ids' must be a character vector with no NAs")
    }

    ## Download and prepare the 'transcripts' data frame.
    attributes <- c("ensembl_transcript_id",
                    "chromosome_name",
                    "strand",
                    "transcript_start",
                    "transcript_end")
    bm_table <- getBM(attributes, filters=filters, values=values, mart=mart)
    transcripts_tx_id <- seq_len(nrow(bm_table))
    transcripts_tx_name <- bm_table$ensembl_transcript_id
    if (any(duplicated(transcripts_tx_name)))
        stop("the 'ensembl_transcript_id' field contains duplicates")
    transcripts <- data.frame(
        tx_id=transcripts_tx_id,
        tx_name=transcripts_tx_name,
        tx_chrom=bm_table$chromosome_name,
        tx_strand=ifelse(bm_table$strand == 1, "+", "-"),
        tx_start=bm_table$transcript_start,
        tx_end=bm_table$transcript_end
    )

    ## Download and prepare the 'splicings' data frame.
    ## Ironically the cds_start and cds_end attributes that we get from
    ## BioMart are pretty useless because they are relative to the coding
    ## mRNA. However, the utr coordinates are relative to the chromosome so
    ## we use them to infer the cds coordinates. We also retrieve the
    ## cds_length attribute as a sanity check.
    attributes <- c("ensembl_transcript_id",
                    "strand",
                    "rank",
                    "ensembl_exon_id",
                    "exon_chrom_start",
                    "exon_chrom_end",
                    "5_utr_start",
                    "5_utr_end",
                    "3_utr_start",
                    "3_utr_end",
                    #"cds_start",
                    #"cds_end",
                    "cds_length")
    bm_table <- getBM(attributes, filters=filters, values=values, mart=mart)
    splicings_tx_id <- as.integer(factor(bm_table$ensembl_transcript_id,
                                         levels=transcripts_tx_name))
    cds_ranges <- .extractCdsRangesFromBiomartTable(bm_table)
    cds_start <- start(cds_ranges)
    cds_start[width(cds_ranges) == 0L] <- NA_integer_
    cds_end <- end(cds_ranges)
    cds_end[width(cds_ranges) == 0L] <- NA_integer_
    splicings <- data.frame(
        tx_id=splicings_tx_id,
        exon_rank=bm_table$rank,
        exon_name=bm_table$ensembl_exon_id,
        exon_start=bm_table$exon_chrom_start,
        exon_end=bm_table$exon_chrom_end,
        cds_start=cds_start,
        cds_end=cds_end
    )

    ## Download and prepare the 'genes' data frame.
    attributes <- c("ensembl_gene_id", "ensembl_transcript_id")
    bm_table <- getBM(attributes, filters=filters, values=values, mart=mart)
    genes_tx_id <- as.integer(factor(bm_table$ensembl_transcript_id,
                                     levels=transcripts_tx_name))
    genes <- data.frame(
        tx_id=genes_tx_id,
        gene_id=bm_table$ensembl_gene_id
    )

    ## Prepare the 'metadata' data frame.
    metadata <- data.frame(
        name=c("Data source",
               "BioMart database",
               "BioMart dataset",
               "BioMart dataset description",
               "BioMart dataset version",
               "Full dataset"),
        value=c("BioMart",
                biomart,
                dataset,
                description,
                version,
                ifelse(is.null(transcript_ids), "yes", "no"))
    )

    ## Call makeTranscriptDb().
    makeTranscriptDb(transcripts, splicings, genes=genes, metadata=metadata)
}

