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

## helper to extract the organism (as Genus and Species) from the dataset
## string.
.extractOrganismFromDatasetDesc <- function(description){
  vals <- unlist(strsplit(description, " "))
  paste(vals[[1]], vals[[2]])
}


.getBiomartDbVersion <- function(biomart, host, port)
{
    marts <- listMarts(mart=biomart, host=host, port=port)

    mart_rowidx <- which(as.character(marts$biomart) == biomart)
    ## This should never happen.
    if (length(mart_rowidx) != 1L)
        stop("found 0 or more than 1 \"", biomart, "\" BioMart database")
    as.character(marts$version)[mart_rowidx]
}

.extractEnsemblReleaseFromDbVersion <- function(db_version)
    sub("^ENSEMBL GENES ([^[:space:]]+) \\(SANGER UK\\)", "\\1", db_version)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'transcripts' data frame.
###

.makeBiomartTranscripts <- function(filters, values, mart, transcript_ids,
                                    biomartAttribGroups,
                                    id_prefix="ensembl_")
{
    message("Download and preprocess the 'transcripts' data frame ... ",
            appendLF=FALSE)
    bm_table <- getBM(biomartAttribGroups[['A1']], filters=filters,
                      values=values, mart=mart, bmHeader=FALSE)
    tx_id_colname <- paste0(id_prefix, "transcript_id")
    if (!is.null(transcript_ids)) {
        idx <- !(transcript_ids %in% bm_table[[tx_id_colname]])
        if (any(idx)) {
            bad_ids <- transcript_ids[idx]
            stop("invalid transcript ids: ",
                 paste0(bad_ids, collapse=", "))
        }
    }
    ## Those are the strictly required fields.
    transcripts0 <- data.frame(
        tx_id=integer(0),
        tx_chrom=character(0),
        tx_strand=character(0),
        tx_start=integer(0),
        tx_end=integer(0)
    )
    if (nrow(bm_table) == 0L) {
        message("OK")
        return(transcripts0)
    }
    tx_id <- seq_len(nrow(bm_table))
    tx_name <- bm_table[[tx_id_colname]]
    ##if (any(duplicated(tx_name)))
    ##    stop(paste("the '",
    ##               tx_id_colname,
    ##               "'transcript_id' attribute contains duplicated values"))
    if (any(duplicated(bm_table)))
      stop("The 'transcripts' data frame from biomart contains duplicated rows.")
    tx_chrom <- as.character(bm_table$chromosome_name)
    tx_strand <- ifelse(bm_table$strand == 1, "+", "-")
    tx_start <- bm_table$transcript_start
    tx_end <- bm_table$transcript_end
    transcripts <- data.frame(
        tx_id=tx_id,
        tx_name=tx_name,
        tx_chrom=tx_chrom,
        tx_strand=tx_strand,
        tx_start=tx_start,
        tx_end=tx_end
    )
    message("OK")
    transcripts
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'chrominfo' data frame.
###

### Returns NULL if it fails to fetch the chromosome lengths from the
### remote resource.
.makeBiomartChrominfo <- function(mart, extra_seqnames=NULL,
                                  circ_seqs=character(0), host, port)
{
    biomart <- biomaRt:::martBM(mart)
    dataset <- biomaRt:::martDataset(mart)
    if (biomart == "ensembl" || substr(biomart, 1, 12) == "plants_mart_") {
        message("Download and preprocess the 'chrominfo' data frame ... ",
                appendLF=FALSE)
        if (biomart == "ensembl") {
            db_version <- .getBiomartDbVersion(biomart, host, port)
            ensembl_release <- .extractEnsemblReleaseFromDbVersion(db_version)
            chromlengths <- try(fetchChromLengthsFromEnsembl(dataset,
                                    release=ensembl_release,
                                    extra_seqnames=extra_seqnames),
                                silent=TRUE)
        } else {
            chromlengths <- try(fetchChromLengthsFromEnsemblPlants(dataset,
                                    extra_seqnames=extra_seqnames),
                                silent=TRUE)
        }
        if (is(chromlengths, "try-error")) {
            message("FAILED! (=> skipped)")
            return(NULL)
        }
        chrominfo <- data.frame(
            chrom=chromlengths$name,
            length=chromlengths$length,
            is_circular=matchCircularity(chromlengths$name, circ_seqs)
        )
        message("OK")
        return(chrominfo)
    }
    NULL
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Allow users to discover 'chrominfo' data frame.
###

getChromInfoFromBiomart <- function(biomart="ensembl",
                                    dataset="hsapiens_gene_ensembl",
                                    id_prefix="ensembl_",
                                    host="www.biomart.org",
                                    port=80)
{
    biomartAttribGroups <- getBiomartAttribGroups(id_prefix)

    mart <- .parseBMMartParams(biomart=biomart,
                               dataset=dataset,
                               host=host,
                               port=port)

    filters <- .parseBMFiltersParams(transcript_ids=NULL, id_prefix)
    values <- .parseBMValuesParams(transcript_ids=NULL)
    transcripts <- .makeBiomartTranscripts(filters, values, mart,
                                           transcript_ids=NULL,
                                           biomartAttribGroups,
                                           id_prefix)
    chrominfo <- .makeBiomartChrominfo(mart,
                                       extra_seqnames=transcripts$tx_chrom,
                                       host=host, port=port)
    chrominfo[,1:2]
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'splicings' data frame.
###

.normUtrCoords <- function(coords)
{
    if (is.numeric(coords))
        return(coords)
    if (is.logical(coords) && all(is.na(coords)))
        return(as.integer(coords))
    stop("BioMart data anomaly: utr coordinates don't ",
         "have a numeric type")
}

.generateBioMartDataAnomalyReport <- function(bm_table, idx, id_prefix, msg)
{
    ## Part 3.
    tx_id_colname <- paste0(id_prefix, "transcript_id")
    tx_ids <- bm_table[[tx_id_colname]]
    first_tx_ids <- unique(tx_ids[idx])
    total_nb_tx <- length(first_tx_ids)
    first_six_only <- total_nb_tx > 6L
    if (first_six_only)
        first_tx_ids <- first_tx_ids[1:6]
    bm_table <- bm_table[tx_ids %in% first_tx_ids, , drop=FALSE]
    bm_table0 <- bm_table[-match(tx_id_colname, names(bm_table))]
    f <- factor(bm_table[[tx_id_colname]], levels=first_tx_ids)
    first_tx_tables <- split(bm_table0, f)
    .DETAILS_INDENT <- "     "
    options(width=getOption("width")-nchar(.DETAILS_INDENT))
    part3 <- lapply(seq_len(length(first_tx_tables)),
                    function(i) {
                        tx_table <- first_tx_tables[[i]]
                        if ("rank" %in% colnames(tx_table)) {
                            oo <- order(tx_table[["rank"]])
                            tx_table <- tx_table[oo, , drop=FALSE]
                        }
                        row.names(tx_table) <- NULL
                        subtitle <- paste0("  ", i, ". Transcript ",
                                           names(first_tx_tables)[i],
                                           ":")
                        details <- capture.output(print(tx_table))
                        c(subtitle, paste0(.DETAILS_INDENT, details))
                    })
    options(width=getOption("width")+nchar(.DETAILS_INDENT))
    part3 <- unlist(part3, use.names=FALSE)
    if (first_six_only)
        part3 <- c(paste("  (Showing only the first 6 out of",
                          total_nb_tx,
                          "transcripts.)"),
                   part3)

    ## Part 1.
    part1 <- "BioMart data anomaly: in the following transcripts, "

    ## Part 2.
    msg[length(msg)] <- paste0(msg[length(msg)], ".")
    part2 <- paste0("  ", msg)

    ## Assemble the parts.
    paste(c(part1, part2, part3), collapse="\n")
}

.stopWithBioMartDataAnomalyReport <- function(bm_table, idx, id_prefix, msg)
{
    msg <- .generateBioMartDataAnomalyReport(bm_table, idx, id_prefix, msg)
    new_length <- nchar(msg) + 5L
    ## 8170L seems to be the maximum possible value for the 'warning.length'
    ## option on my machine (R-2.15 r58124, 64-bit Ubuntu).
    if (new_length > 8170L)
        new_length <- 8170L
    if (new_length >= getOption("warning.length")) {
        old_length <- getOption("warning.length")
        on.exit(options(warning.length=old_length))
        options(warning.length=new_length)
    }
    stop(msg)
}

.warningWithBioMartDataAnomalyReport <- function(bm_table, idx, id_prefix, msg)
{
    msg <- .generateBioMartDataAnomalyReport(bm_table, idx, id_prefix, msg)
    new_length <- nchar(msg) + 5L
    ## 8170L seems to be the maximum possible value for the 'warning.length'
    ## option on my machine (R-2.15 r58124, 64-bit Ubuntu).
    if (new_length > 8170L)
        new_length <- 8170L
    if (new_length >= getOption("warning.length")) {
        old_length <- getOption("warning.length")
        on.exit(options(warning.length=old_length))
        options(warning.length=new_length)
    }
    warning(msg)
}

.utrIsNa <- function(utr_start, utr_end, exon_start, exon_end,
                     what_utr, bm_table, id_prefix="ensembl_")
{
    is_na <- is.na(utr_start)
    if (!identical(is_na, is.na(utr_end)))
        stop("BioMart data anomaly: ",
             "NAs in \"", what_utr, "_utr_start\" attribute don't match ",
             "NAs in \"", what_utr, "_utr_end\" attribute")
    idx <- which(utr_start > utr_end + 1L)
    if (length(idx) != 0L) {
        msg <- paste0("the ", what_utr, "' UTRs have a start > end")
        .stopWithBioMartDataAnomalyReport(bm_table, idx, id_prefix, msg)
    }
    idx <- which(utr_start < exon_start | exon_end < utr_end)
    if (length(idx) != 0L) {
        msg <- paste0("the ", what_utr, "' UTRs ",
                      "are not within the exon limits")
        .stopWithBioMartDataAnomalyReport(bm_table, idx, id_prefix, msg)
    }
    is_na | utr_start == utr_end + 1L
}

.extractCdsRangesFromBiomartTable <- function(bm_table, id_prefix="ensembl_")
{
    if (nrow(bm_table) == 0L)
        return(IRanges())
    strand <- bm_table[["strand"]]
    cds_start <- exon_start <- bm_table[["exon_chrom_start"]]
    cds_end <- exon_end <- bm_table[["exon_chrom_end"]]
    utr5_start <- .normUtrCoords(bm_table[["5_utr_start"]])
    utr5_end <- .normUtrCoords(bm_table[["5_utr_end"]])
    utr3_start <- .normUtrCoords(bm_table[["3_utr_start"]])
    utr3_end <- .normUtrCoords(bm_table[["3_utr_end"]])

    if (!all(strand %in% c(1, -1)))
        stop("BioMart data anomaly: \"strand\" attribute should be 1 or -1")
    if (!is.numeric(exon_start) || !is.numeric(exon_end))
        stop("BioMart data anomaly: exon coordinates don't ",
             "have a numeric type")

    no_utr5 <- .utrIsNa(utr5_start, utr5_end, exon_start, exon_end,
                        "5", bm_table, id_prefix)
    no_utr3 <- .utrIsNa(utr3_start, utr3_end, exon_start, exon_end,
                        "3", bm_table, id_prefix)

    idx <- strand == 1 & !no_utr5
    if (!all(utr5_start[idx] == exon_start[idx])) {
        msg <- c("located on the plus strand, the 5' UTRs don't start",
                 "where their corresponding exon starts")
        .stopWithBioMartDataAnomalyReport(bm_table, idx, id_prefix, msg)
    }
    cds_start[idx] <- utr5_end[idx] + 1L

    idx <- strand == 1 & !no_utr3
    if (!all(utr3_end[idx] == exon_end[idx])) {
        msg <- c("located on the plus strand, the 3' UTRs don't end",
                 "where their corresponding exon ends")
        .stopWithBioMartDataAnomalyReport(bm_table, idx, id_prefix, msg)
    }
    cds_end[idx] <- utr3_start[idx] - 1L

    idx <- strand == -1 & !no_utr3
    if (!all(utr3_start[idx] == exon_start[idx])) {
        msg <- c("located on the minus strand, the 3' UTRs don't start",
                 "where their corresponding exon starts")
        .stopWithBioMartDataAnomalyReport(bm_table, idx, id_prefix, msg)
    }
    cds_start[idx] <- utr3_end[idx] + 1L

    idx <- strand == -1 & !no_utr5
    if (!all(utr5_end[idx] == exon_end[idx])) {
        msg <- c("located on the minus strand, the 5' UTRs don't end",
                 "where their corresponding exon ends")
        .stopWithBioMartDataAnomalyReport(bm_table, idx, id_prefix, msg)
    }
    cds_end[idx] <- utr5_start[idx] - 1L

    ans <- IRanges(start=cds_start, end=cds_end)
    if (length(ans) != 0L) {
        tx_id_colname <- paste0(id_prefix, "transcript_id")
        cds_cumlength <-
            sapply(split(width(ans), bm_table[[tx_id_colname]]), sum)
        idx <- which(cds_cumlength[as.vector(bm_table[[tx_id_colname]])] !=
                     bm_table$cds_length)
        if (length(idx) != 0L) {
            msg <- c("the CDS total length inferred from the exon and UTR info",
                     "doesn't match the \"cds_length\" attribute from BioMart")
            .warningWithBioMartDataAnomalyReport(bm_table, idx, id_prefix, msg)
        }
        #idx <- which(cds_cumlength %% 3L != 0L)
        #if (length(idx) != 0L) {
        #    msg <- c("the CDS total length (\"cds_length\" attribute) is not",
        #             "a multiple of 3")
        #    .warningWithBioMartDataAnomalyReport(bm_table, idx, id_prefix, msg)
        #}
    }
    ans
}

.makeCdsDataFrameFromRanges <- function(cds_ranges)
{
    nocds_idx <- width(cds_ranges) == 0L
    cds_start <- start(cds_ranges)
    cds_start[nocds_idx] <- NA_integer_
    cds_end <- end(cds_ranges)
    cds_end[nocds_idx] <- NA_integer_
    data.frame(cds_start=cds_start, cds_end=cds_end)
}

### Surprisingly the cds_start and cds_end attributes that we get from
### BioMart are pretty useless because they are relative to the coding
### mRNA. However, the utr coordinates are relative to the chromosome so
### we use them to infer the cds coordinates. We also retrieve the
### cds_length attribute to do a sanity check.
.makeBiomartSplicings <- function(filters, values, mart, transcripts_tx_id,
                                  biomartAttribGroups, id_prefix="ensembl_")
{

    ## Those are the strictly required fields.
    splicings0 <- data.frame(
        tx_id=integer(0),
        exon_rank=integer(0),
        exon_start=integer(0),
        exon_end=integer(0)
    )
    if (length(transcripts_tx_id) == 0L)
        return(splicings0)
    message("Download and preprocess the 'splicings' data frame ... ",
            appendLF=FALSE)
    allattribs <- listAttributes(mart)$name
    attributes <- biomartAttribGroups[['A2']]
    if (all(biomartAttribGroups[['B']] %in% allattribs))
      attributes <- c(attributes, biomartAttribGroups[['B']])
    if (all(biomartAttribGroups[['C1']] %in% allattribs))
      attributes <- c(attributes, biomartAttribGroups[['C1']])
    if (all(biomartAttribGroups[['C2']] %in% allattribs))
      attributes <- c(attributes, biomartAttribGroups[['C2']])
    if ("cds_length" %in% allattribs)
        attributes <- c(attributes, "cds_length")
    bm_table <- getBM(attributes, filters=filters, values=values, mart=mart,
                      bmHeader=FALSE)
    tx_id_colname <- paste0(id_prefix, "transcript_id")
    splicings_tx_id <- transcripts_tx_id[bm_table[[tx_id_colname]]]
    exon_id_col_name <- paste0(id_prefix, "exon_id")
    splicings <- data.frame(
        tx_id=splicings_tx_id,
        exon_rank=bm_table$rank,
        exon_name=bm_table[[exon_id_col_name]],
        exon_start=bm_table$exon_chrom_start,
        exon_end=bm_table$exon_chrom_end
    )
    if ((all(biomartAttribGroups[['C1']] %in% allattribs) ||
         all(biomartAttribGroups[['C2']] %in% allattribs))
     && ("cds_length" %in% allattribs)) {
        cds_ranges <- .extractCdsRangesFromBiomartTable(bm_table, id_prefix)
        splicings <- cbind(splicings, .makeCdsDataFrameFromRanges(cds_ranges))
    }
    message("OK")
    splicings  
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'genes' data frame.
###

.makeBiomartGenes <- function(filters, values, mart,
                              transcripts_tx_id, biomartAttribGroups,
                              id_prefix="ensembl_")
{
    message("Download and preprocess the 'genes' data frame ... ",
            appendLF=FALSE)
    attributes <- c(biomartAttribGroups[['G']],
                    paste0(id_prefix, "transcript_id"))
    bm_table <- getBM(attributes, filters=filters, values=values, mart=mart,
                      bmHeader=FALSE)
    tx_id_colname <- paste0(id_prefix, "transcript_id")
    genes_tx_id <- transcripts_tx_id[bm_table[[tx_id_colname]]]
    message("OK")
    gene_id_col_name <- paste0(id_prefix, "gene_id")
    data.frame(
        tx_id=genes_tx_id,
        gene_id=bm_table[[gene_id_col_name]]
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Prepare the 'metadata' data frame.
###

.prepareBiomartMetadata <- function(mart, is_full_dataset, host, port,
                                    miRBaseBuild)
{
    message("Prepare the 'metadata' data frame ... ",
            appendLF=FALSE)
    biomart <- biomaRt:::martBM(mart)
    dataset <- biomaRt:::martDataset(mart)
    mart_url <- biomaRt:::martHost(mart)
    mart_url <- sub("^[^/]+//", "", mart_url)
    mart_url <- unlist(strsplit(mart_url, "/"))[1]
    db_version <- .getBiomartDbVersion(biomart, host, port)
    datasets <- listDatasets(mart)
    dataset_rowidx <- which(as.character(datasets$dataset) == dataset)
    ## This should never happen (the above call to useMart() would have failed
    ## in the first place).
    if (length(dataset_rowidx) != 1L)
        stop("the BioMart database \"", biomaRt:::martBM(mart),
             "\" has no (or more than one) \"", dataset, "\" datasets")
    description <- as.character(datasets$description)[dataset_rowidx]
    dataset_version <- as.character(datasets$version)[dataset_rowidx]
    organism <- .extractOrganismFromDatasetDesc(description)
    if (!isSingleStringOrNA(miRBaseBuild))
        stop("'miRBaseBuild' must be a a single string or NA")
    message("OK")
    data.frame(
        name=c("Data source",
               "Organism",
               "Resource URL",
               "BioMart database",
               "BioMart database version",
               "BioMart dataset",
               "BioMart dataset description",
               "BioMart dataset version",
               "Full dataset",
               "miRBase build ID"),
        value=c("BioMart",
                organism,
                mart_url,
                biomart,
                db_version,
                dataset,
                description,
                dataset_version,
                ifelse(is_full_dataset, "yes", "no"),
                miRBaseBuild)
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTranscriptDbFromBiomart()
###

.parseBMMartParams <- function(biomart="ensembl",
                               dataset="hsapiens_gene_ensembl",
                               host, port)
{
    if (is.factor(biomart))
        biomart <- as.character(biomart)
    if (is(dataset, "AsIs"))
        dataset <- as.character(dataset)
    if (!isSingleString(biomart))
        stop("'biomart' must be a single string")
    useMart(biomart=biomart, dataset=dataset, host=host, port=port)
}

.parseBMFiltersParams <- function(transcript_ids, id_prefix="ensembl_")
{
    if (is.null(transcript_ids)) {
        filters <- ""
    } else if (is.character(transcript_ids)
            && !any(is.na(transcript_ids))) {
        filters <- paste0(id_prefix, "transcript_id")
    }
    filters
}

.parseBMValuesParams <- function(transcript_ids)
{
    if (is.null(transcript_ids)) {
        values <- ""
    }else if (is.character(transcript_ids)
            && !any(is.na(transcript_ids))) {
        if (length(transcript_ids) == 0L)
            values <- "____a_very_unlikely_valid_transcript_id____"
        else
            values <- transcript_ids
    } else {
        stop("'transcript_ids' must be a character vector with no NAs")
    }
    values
}


## .testMakeTxDbFromBMParams <- function(biomart="ensembl",
##                                       dataset="hsapiens_gene_ensembl",
##                                       circ_seqs=DEFAULT_CIRC_SEQS,
##                                       transcript_ids=NULL)
## {
    ## if (is.factor(biomart))
    ##     biomart <- as.character(biomart)
    ## if (is(dataset, "AsIs"))
    ##     dataset <- as.character(dataset)
    ## if (!isSingleString(biomart))
    ##     stop("'biomart' must be a single string")
    ## mart <- useMart(biomart=biomart, dataset=dataset)
  
    ## if (is.null(transcript_ids)) {
    ##     filters <- values <- ""
    ## } else if (is.character(transcript_ids)
    ##         && !any(is.na(transcript_ids))) {
    ##     filters <- "ensembl_transcript_id" 
    ##     if (length(transcript_ids) == 0L)
    ##         values <- "____a_very_unlikely_valid_transcript_id____"
    ##     else
    ##         values <- transcript_ids
    ## } else {
    ##     stop("'transcript_ids' must be a character vector with no NAs")
    ## }
## }


### Note that listMarts() and listDatasets() are returning data frames where
### the columns are character factors for the former and "AsIs" character
### vectors for the latter.
makeTranscriptDbFromBiomart <- function(biomart="ensembl",
                                        dataset="hsapiens_gene_ensembl",
                                        transcript_ids=NULL,
                                        circ_seqs=DEFAULT_CIRC_SEQS,
                                        filters="",
                                        id_prefix="ensembl_",
                                        host="www.biomart.org",
                                        port=80,
                                        miRBaseBuild=NA)
{
    ## Could be that the user got the 'biomart' and/or 'dataset' values
    ## programmatically via calls to listMarts() and/or listDatasets().
    mart <- .parseBMMartParams(biomart=biomart,
                               dataset=dataset,
                               host=host,
                               port=port)

    ##combines user-specified filters with those created from supplied
    ##transcripts_ids
    .mergeTxIDsAndFilters <- function(transcript_ids, filters, id_prefix) {
      if (filters == "")
        filters <- list()

      if (class(filters) != "list")
        stop("filters parameter must be a named list")

      if(length(filters) != 0) {
        if(is.null(names(filters)))
          stop("filters parameter must be a named list")
      }
      
      transcript_filters <- .parseBMFiltersParams(transcript_ids, id_prefix)
      transcript_values <- .parseBMValuesParams(transcript_ids)

      ##merge transcript_ids into filters
      transcript_list <- list()
      if(transcript_filters != "") {
        transcript_list <- list(transcript_values)
        names(transcript_list) <- transcript_filters
      }

      transcripts_and_filters <- append(filters, transcript_list)

      f <- ""
      v <- ""
      if (length(transcripts_and_filters) > 0) {
        f <- names(transcripts_and_filters)
        v <- unname(transcripts_and_filters)
      }
      res <- list()
      res[['filters']] <- f
      res[['values']] <- v
      return(res)
    }
    
    return_list <- .mergeTxIDsAndFilters(transcript_ids,
                                         filters, id_prefix)
    filters <- return_list$filters
    values <- return_list$values
    
    biomartAttribGroups <- getBiomartAttribGroups(id_prefix)
    transcripts <- .makeBiomartTranscripts(filters, values, mart,
                                           transcript_ids,
                                           biomartAttribGroups,
                                           id_prefix)
    transcripts_tx_id <- transcripts$tx_id
    names(transcripts_tx_id) <- transcripts$tx_name
    chrominfo <- .makeBiomartChrominfo(mart,
                                       extra_seqnames=transcripts$tx_chrom,
                                       circ_seqs=circ_seqs,
                                       host, port)
    splicings <- .makeBiomartSplicings(filters, values, mart,
                                       transcripts_tx_id,
                                       biomartAttribGroups,
                                       id_prefix=id_prefix)
    genes <- .makeBiomartGenes(filters, values, mart, transcripts_tx_id,
                               biomartAttribGroups, id_prefix)
    metadata <- .prepareBiomartMetadata(mart, is.null(transcript_ids), host,
                                        port, miRBaseBuild)

    message("Make the TxDb object ... ", appendLF=FALSE)
    txdb <- makeTranscriptDb(transcripts, splicings,
                             genes=genes, chrominfo=chrominfo,
                             metadata=metadata, reassign.ids=TRUE)
    message("OK")
    txdb
}

