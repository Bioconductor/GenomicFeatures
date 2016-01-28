### =========================================================================
### makeTxDbFromBiomart()
### -------------------------------------------------------------------------
###
### For people who want to tap BioMart.
### Typical use:
###   txdb <- makeTxDbFromBiomart("hsapiens_gene_ensembl")
### Speed:
###   - for biomart="ENSEMBL_MART_ENSEMBL" and dataset="hsapiens_gene_ensembl":
###       (1) download takes about 8 min.
###       (2) db creation takes about 60-65 sec.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some helper functions to facilitate working with the biomaRt package.
###

### A thin wrapper to useMart() that checks the user-supplied arguments.
.useMart2 <- function(biomart="ENSEMBL_MART_ENSEMBL",
                      dataset="hsapiens_gene_ensembl",
                      host="www.ensembl.org",
                      port=80)
{
    ### Could be that the user got the 'biomart' and/or 'dataset' values
    ### programmatically via calls to listMarts() and/or listDatasets(). Note
    ### that listMarts() and listDatasets() are returning data frames where
    ### the columns are factors for the former and "AsIs" character vectors
    ### for the latter.
    if (is.factor(biomart))
        biomart <- as.character(biomart)
    if (!(isSingleString(biomart) && biomart != ""))
        stop("'biomart' must be a single non-empty string")
    if (is(dataset, "AsIs"))
        dataset <- as.character(dataset)
    if (!(isSingleString(dataset) && dataset != ""))
        stop("'dataset' must be a single non-empty string")
    if (!(isSingleString(host) && host != ""))
        stop("'host' must be a single non-empty string")
    if (!(isSingleNumber(port) && port > 0))
        stop("'port' must be a single positive number")
    useMart(biomart=biomart, dataset=dataset, host=host, port=port)
}

### TODO: Share this with normalization of 'vals' arg in the transcripts(),
### exons(), cds(), and genes() extractors.
.normarg_filters <- function(filters)
{
    if (is.null(filters) || identical(filters, ""))
        return(setNames(list(), character(0)))
    if (!is.list(filters))
        stop("'filters' must be a named list")
    if (length(filters) == 0L)
        return(setNames(list(), character(0)))
    filters_names <- names(filters)
    if (is.null(filters_names))
        stop("'filters' must be a named list")
    if (any(filters_names %in% c(NA, "")))
        stop("names on 'filters' cannot be NA or the empty string")
    if (anyDuplicated(filters_names))
        stop("names on 'filters' must be unique")
    if (!all(sapply(filters, is.atomic)))
        stop("'filters' list elements must be atomic")
    if (any(sapply(filters, anyNA)))
        stop("'filters' list elements cannot contain NAs")
    filters
}

### A thin wrapper to getBM() that takes the filters in the form of a named
### list.
.getBM2 <- function(attributes, filters=NULL, ...)
{
    filters <- .normarg_filters(filters)
    if (length(filters) == 0L) {
        bm_filters <- bm_values <- ""
    } else {
        bm_filters <- names(filters)
        bm_values <- unname(filters)
        bm_values[elementNROWS(bm_values) == 0L] <- paste0(
            "____this_is_a_very_unlikely_valid_value_but_you_never_know_",
            "this_is_just_a_dirty_hack_to_work_around_getBM_",
            "misinterpretation_of_empty_list_elements_in_values____")
    }
    getBM(attributes, filters=bm_filters, values=bm_values, ...)
}

.normarg_id_prefix <- function(id_prefix)
{
    if (!isSingleString(id_prefix))
        stop("'id_prefix' must be a single string")
    id_prefix
}

### Add filter created from user-supplied transcript_ids to user-specified
### filters.
.add_tx_id_filter <- function(filters,
                              transcript_ids=NULL, id_prefix="ensembl_")
{
    filters <- .normarg_filters(filters)
    tx_name_colname <- paste0(id_prefix, "transcript_id")
    if (is.null(transcript_ids)) {
        if (tx_name_colname %in% names(filters))
            warning(wmsg("transcript ids should be specified via the ",
                         "'transcript_ids' rather than the 'filters' argument"))
            return(filters)
    }
    if (!is.character(transcript_ids))
        stop("'transcript_ids' must ba a character vector")
    if (any(is.na(transcript_ids)))
        stop("'transcript_ids' cannot contain NAs")
    if (tx_name_colname %in% names(filters))
        stop(wmsg("transcript ids cannot be specified via the ",
                  "'transcript_ids' and 'filters' arguments ",
                  "at the same time"))
    filters[[tx_name_colname]] <- transcript_ids
    filters
}

## helper to extract the organism (as Genus and Species) from the dataset
## string.
.extractOrganismFromDatasetDesc <- function(description){
  vals <- unlist(strsplit(description, " "))
  paste(vals[[1]], vals[[2]])
}

.getBiomartDbVersion <- function(mart, host, port, biomart) 
{
    marts <- listMarts(mart=mart, host=host, port=port)

    mart_rowidx <- which(as.character(marts$biomart) == biomart)
    ## This should never happen.
    if (length(mart_rowidx) != 1L)
        stop("found 0 or more than 1 \"", biomart, "\" BioMart database")
    as.character(marts$version)[mart_rowidx]
}

.extractEnsemblReleaseFromDbVersion <- function(db_version)
{
    db_version <- tolower(db_version)
    #sub("^ensembl genes ([^[:space:]]+) \\(sanger uk\\)", "\\1", db_version)
    sub("^ensembl genes ([0-9]+).*$", "\\1", db_version)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'transcripts' data frame.
###

.makeBiomartTranscripts <- function(filters, mart, transcript_ids,
                                    recognized_attribs,
                                    id_prefix="ensembl_")
{
    message("Download and preprocess the 'transcripts' data frame ... ",
            appendLF=FALSE)
    bm_result <- .getBM2(recognized_attribs[['T']], filters,
                         mart=mart, bmHeader=FALSE)
    tx_name_colname <- paste0(id_prefix, "transcript_id")
    tx_name <- bm_result[[tx_name_colname]]
    if (!is.null(transcript_ids)) {
        idx <- !(transcript_ids %in% tx_name)
        if (any(idx)) {
            bad_ids <- transcript_ids[idx]
            stop(wmsg("invalid transcript ids: ",
                      paste0(bad_ids, collapse=", ")))
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
    if (nrow(bm_result) == 0L) {
        message("OK")
        return(transcripts0)
    }
    tx_id <- seq_len(nrow(bm_result))
    ##if (any(duplicated(tx_name)))
    ##    stop(wmsg("the '",
    ##              tx_name_colname,
    ##              "'transcript_id' attribute contains duplicated values"))
    if (any(duplicated(bm_result)))
        stop(wmsg("the 'transcripts' data frame obtained from biomart ",
                  "contains duplicated rows"))
    tx_type <- as.character(bm_result$transcript_biotype)
    tx_chrom <- as.character(bm_result$chromosome_name)
    tx_strand <- ifelse(bm_result$strand == 1, "+", "-")
    tx_start <- bm_result$transcript_start
    tx_end <- bm_result$transcript_end
    transcripts <- data.frame(
        tx_id=tx_id,
        tx_name=tx_name,
        tx_type=tx_type,
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
                                  circ_seqs=DEFAULT_CIRC_SEQS, host, port)
{
    biomart <- biomaRt:::martBM(mart)
    dataset <- biomaRt:::martDataset(mart)
    is_ensembl_mart <- tolower(substr(biomart, 1, 7)) == "ensembl"
    is_plants_mart <- tolower(substr(biomart, 1, 11)) == "plants_mart"
    if (is_ensembl_mart || is_plants_mart) {
        message("Download and preprocess the 'chrominfo' data frame ... ",
                appendLF=FALSE)
        if (is_ensembl_mart) {
            if (tolower(host) == "grch37.ensembl.org") {
                ## Ensembl GRCh37 mart
                chromlengths <- try(fetchChromLengthsFromEnsembl(dataset,
                                        use.grch37=TRUE,
                                        extra_seqnames=extra_seqnames),
                                    silent=TRUE)
            } else {
                ## Ensembl mart
                db_version <- .getBiomartDbVersion(mart, host, port, biomart)
                ensembl_release <-
                              .extractEnsemblReleaseFromDbVersion(db_version)
                chromlengths <- try(fetchChromLengthsFromEnsembl(dataset,
                                        release=ensembl_release,
                                        extra_seqnames=extra_seqnames),
                                    silent=TRUE)
            }
        } else {
            ## Plants mart
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
            is_circular=make_circ_flags_from_circ_seqs(chromlengths$name,
                                                       circ_seqs)
        )
        message("OK")
        return(chrominfo)
    }
    NULL
}

## User-friendly wrapper to .makeBiomartChrominfo().
getChromInfoFromBiomart <- function(biomart="ENSEMBL_MART_ENSEMBL",
                                    dataset="hsapiens_gene_ensembl",
                                    id_prefix="ensembl_",
                                    host="www.ensembl.org",
                                    port=80)
{
    mart <- .useMart2(biomart=biomart, dataset=dataset, host=host, port=port)
    id_prefix <- .normarg_id_prefix(id_prefix)
    recognized_attribs <- recognizedBiomartAttribs(id_prefix)
    transcripts <- .makeBiomartTranscripts(NULL, mart,
                                           transcript_ids=NULL,
                                           recognized_attribs,
                                           id_prefix)
    chrominfo <- .makeBiomartChrominfo(mart,
                                       extra_seqnames=transcripts$tx_chrom,
                                       host=host, port=port)
    chrominfo[ , 1:2, drop=FALSE]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'splicings' data frame.
###

.extract_numeric_attrib <- function(bm_result, attrib)
{
    ans <- bm_result[[attrib]]
    if (is.numeric(ans))
        return(ans)
    if (is.logical(ans) && all(is.na(ans)))
        return(as.integer(ans))
    stop(wmsg("BioMart fatal data anomaly: ",
              "\"", attrib, "\" attribute is not numeric"))
}

.generate_BioMart_data_anomaly_report <- function(error_type, bm_result, idx,
                                                  id_prefix, msg)
{
    ## Part 3.
    tx_name_colname <- paste0(id_prefix, "transcript_id")
    tx_name <- bm_result[[tx_name_colname]]
    first_tx_names <- unique(tx_name[idx])
    total_nb_tx <- length(first_tx_names)
    first_three_only <- total_nb_tx > 3L
    if (first_three_only)
        first_tx_names <- first_tx_names[1:3]
    bm_result <- S4Vectors:::extract_data_frame_rows(bm_result,
                                     tx_name %in% first_tx_names)
    bm_result0 <- bm_result[-match(tx_name_colname, names(bm_result))]
    f <- factor(bm_result[[tx_name_colname]], levels=first_tx_names)
    first_tx_tables <- split(bm_result0, f)
    .DETAILS_INDENT <- "     "
    options(width=getOption("width")-nchar(.DETAILS_INDENT))
    part3 <- lapply(seq_len(length(first_tx_tables)),
                    function(i) {
                        tx_table <- first_tx_tables[[i]]
                        if ("rank" %in% colnames(tx_table)) {
                            oo <- order(tx_table[["rank"]])
                            tx_table <-
                              S4Vectors:::extract_data_frame_rows(tx_table, oo)
                        } else {
                            rownames(tx_table) <- NULL
                        }
                        subtitle <- paste0("  ", i, ". Transcript ",
                                           names(first_tx_tables)[i],
                                           ":")
                        details <- capture.output(print(tx_table))
                        c(subtitle, paste0(.DETAILS_INDENT, details))
                    })
    options(width=getOption("width")+nchar(.DETAILS_INDENT))
    part3 <- unlist(part3, use.names=FALSE)
    if (first_three_only)
        part3 <- c(paste("  (Showing only the first 3 out of",
                          total_nb_tx,
                          "transcripts.)"),
                   part3)

    ## Part 1.
    part1 <- paste0(error_type, ": in the following transcripts, ")

    ## Part 2.
    msg[length(msg)] <- paste0(msg[length(msg)], ".")
    part2 <- paste0("  ", msg)

    ## Assemble the parts.
    paste(c(part1, part2, part3), collapse="\n")
}

.stop_on_BioMart_data_anomaly <- function(bm_result, idx, id_prefix, msg)
{
    msg <- .generate_BioMart_data_anomaly_report("BioMart fatal data anomaly",
                                                 bm_result, idx, id_prefix, msg)
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

.warning_on_BioMart_data_anomaly <- function(bm_result, idx, id_prefix, msg)
{
    msg <- .generate_BioMart_data_anomaly_report("BioMart data anomaly",
                                                 bm_result, idx, id_prefix, msg)
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

.has_utr <- function(utr_start, utr_end, exon_start, exon_end,
                     what_utr, bm_result, id_prefix="ensembl_")
{
    is_na <- is.na(utr_start)
    if (!identical(is_na, is.na(utr_end)))
        stop(wmsg("BioMart fatal data anomaly: ",
                  "NAs in \"", what_utr, "_utr_start\" attribute don't match ",
                  "NAs in \"", what_utr, "_utr_end\" attribute"))
    idx <- which(utr_start > utr_end + 1L)
    if (length(idx) != 0L) {
        msg <- paste0("the ", what_utr, "' UTRs have a start > end + 1")
        .stop_on_BioMart_data_anomaly(bm_result, idx, id_prefix, msg)
    }
    idx <- which(utr_start < exon_start | exon_end < utr_end)
    if (length(idx) != 0L) {
        msg <- paste0("the ", what_utr, "' UTRs ",
                      "are not within the exon limits")
        .stop_on_BioMart_data_anomaly(bm_result, idx, id_prefix, msg)
    }
    !(is_na | utr_start == utr_end + 1L)
}

.extract_cds_ranges_from_C1 <- function(bm_result, id_prefix="ensembl_")
{
    cds_start <- .extract_numeric_attrib(bm_result, "genomic_coding_start")
    cds_end <- .extract_numeric_attrib(bm_result, "genomic_coding_end")
    is_na <- is.na(cds_start)
    if (!identical(is_na, is.na(cds_end)))
        stop(wmsg("BioMart fatal data anomaly: ",
                  "NAs in \"genomic_coding_start\" attribute don't match ",
                  "NAs in \"genomic_coding_end\" attribute"))

    ## Exons with no CDS get a CDS of width 0.
    no_cds_idx <- which(is_na)
    exon_start <- bm_result[["exon_chrom_start"]]
    cds_start[no_cds_idx] <- exon_start[no_cds_idx]
    cds_end[no_cds_idx] <- cds_start[no_cds_idx] - 1L

    IRanges(start=cds_start, end=cds_end)
}

### These errors in UTR representation are non fatal but trigger rejection of
### the corresponding transcripts with a warning.
.BIOMART_UTR_ERROR <- c(
    "located on the + strand, \"5_utr_start\" must match \"exon_chrom_start\"",
    "located on the + strand, \"3_utr_end\" must match \"exon_chrom_end\"",
    "located on the - strand, \"3_utr_start\" must match \"exon_chrom_start\"",
    "located on the - strand, \"5_utr_end\" must match \"exon_chrom_end\""
)

.warning_on_BioMart_utr_anomaly <- function(bm_result, idx, id_prefix,
                                            utr_anomaly)
{
    msg <- c(.BIOMART_UTR_ERROR[[utr_anomaly]],
             " (these transcripts were dropped)")
    .warning_on_BioMart_data_anomaly(bm_result, idx, id_prefix, msg)
}

.extract_cds_ranges_from_C2 <- function(bm_result, id_prefix="ensembl_")
{
    strand <- bm_result[["strand"]]
    if (!all(strand %in% c(1, -1)))
        stop(wmsg("BioMart fatal data anomaly: ",
                  "\"strand\" attribute should be 1 or -1"))

    cds_start <- exon_start <- bm_result[["exon_chrom_start"]]
    cds_end <- exon_end <- bm_result[["exon_chrom_end"]]
    utr_anomaly <- integer(nrow(bm_result))

    utr5_start <- .extract_numeric_attrib(bm_result, "5_utr_start")
    utr5_end <- .extract_numeric_attrib(bm_result, "5_utr_end")
    utr3_start <- .extract_numeric_attrib(bm_result, "3_utr_start")
    utr3_end <- .extract_numeric_attrib(bm_result, "3_utr_end")

    has_utr5 <- .has_utr(utr5_start, utr5_end, exon_start, exon_end,
                         "5", bm_result, id_prefix)
    has_utr3 <- .has_utr(utr3_start, utr3_end, exon_start, exon_end,
                         "3", bm_result, id_prefix)

    idx <- which(strand == 1 & has_utr5)
    bad_idx <- idx[utr5_start[idx] != exon_start[idx]]
    if (length(bad_idx) != 0L)
        .warning_on_BioMart_utr_anomaly(bm_result, bad_idx, id_prefix,
                                        utr_anomaly[bad_idx] <- 1L)
    cds_start[idx] <- utr5_end[idx] + 1L

    idx <- which(strand == 1 & has_utr3)
    bad_idx <- idx[utr3_end[idx] != exon_end[idx]]
    if (length(bad_idx) != 0L)
        .warning_on_BioMart_utr_anomaly(bm_result, bad_idx, id_prefix,
                                        utr_anomaly[bad_idx] <- 2L)
    cds_end[idx] <- utr3_start[idx] - 1L

    idx <- which(strand == -1 & has_utr3)
    bad_idx <- idx[utr3_start[idx] != exon_start[idx]]
    if (length(bad_idx) != 0L)
        .warning_on_BioMart_utr_anomaly(bm_result, bad_idx, id_prefix,
                                        utr_anomaly[bad_idx] <- 3L)
    cds_start[idx] <- utr3_end[idx] + 1L

    idx <- which(strand == -1 & has_utr5)
    bad_idx <- idx[utr5_end[idx] != exon_end[idx]]
    if (length(bad_idx) != 0L)
        .warning_on_BioMart_utr_anomaly(bm_result, bad_idx, id_prefix,
                                        utr_anomaly[bad_idx] <- 4L)
    cds_end[idx] <- utr5_start[idx] - 1L

    ## Exons with no CDS get a CDS of width 0.
    cds_relative_start <- bm_result[["cds_start"]]
    no_cds_idx <- which(is.na(cds_relative_start))
    cds_end[no_cds_idx] <- cds_start[no_cds_idx] - 1L

    ans <- IRanges(start=cds_start, end=cds_end)
    mcols(ans) <- DataFrame(utr_anomaly=utr_anomaly)
    ans
}

.check_cds <- function(cds_ranges, cds_width, bm_result, id_prefix="ensembl_")
{
    idx <- which(width(cds_ranges) != cds_width)
    if (length(idx) != 0L) {
        msg <- c("the CDS/UTR genomic coordinates are inconsistent with the ",
                 "\"cds_start\" and \"cds_end\" attributes")
        .warning_on_BioMart_data_anomaly(bm_result, idx, id_prefix, msg)
    }

    tx_name_colname <- paste0(id_prefix, "transcript_id")
    tx_name <- bm_result[ , tx_name_colname]
    cds_length2 <- sapply(split(width(cds_ranges), tx_name), sum)
    cds_length2 <- cds_length2[as.character(tx_name)]

    cds_length <- bm_result$cds_length
    if (!is.null(cds_length)) { 
        idx <- which(cds_length2 != cds_length)
        if (length(idx) != 0L) {
            msg <- c("the CDS length inferred from the CDS/UTR genomic ",
                     "coordinates doesn't match the \"cds_length\" attribute")
            .warning_on_BioMart_data_anomaly(bm_result, idx, id_prefix, msg)
        }
    }

    ## Too many transcripts in the ensembl/hsapiens_gene_ensembl dataset don't
    ## pass the sanity check below (20256 transcripts in Ensembl release 75).
    ## This makes makeTxDbFromBiomart() a little bit too noisy so we
    ## comment this out for now.
    #idx <- which(cds_length2 %% 3L != 0L)
    #if (length(idx) != 0L) {
    #    msg <- c("the CDS length inferred from the CDS/UTR genomic ",
    #             "coordinates is not a multiple of 3")
    #    .warning_on_BioMart_data_anomaly(bm_result, idx, id_prefix, msg)
    #}
}

.extract_cds_ranges_from_bm_result <- function(bm_result, id_prefix="ensembl_")
{
    if (nrow(bm_result) == 0L)
        return(IRanges())

    exon_start <- bm_result[["exon_chrom_start"]]
    exon_end <- bm_result[["exon_chrom_end"]]
    if (!is.numeric(exon_start) || !is.numeric(exon_end))
        stop("BioMart data anomaly: \"exon_chrom_start\" and/or ",
             "\"exon_chrom_end\" attributes are not numeric")

    ## BE AWARE that the "cds_start" and "cds_end" attributes that we get
    ## from BioMart are the CDS coordinates relative to the coding mRNA!
    ## See IMPORTANT NOTE ABOUT GROUP D1 in findCompatibleMarts.R for more
    ## information.
    cds_relative_start <- .extract_numeric_attrib(bm_result, "cds_start")
    cds_relative_end <- .extract_numeric_attrib(bm_result, "cds_end")
    is_na <- is.na(cds_relative_start)
    if (!identical(is_na, is.na(cds_relative_end)))
        stop("BioMart data anomaly: ",
             "NAs in \"cds_start\" attribute don't match ",
             "NAs in \"cds_end\" attribute")
    no_cds_idx <- which(is_na)
    cds_width <- cds_relative_end - cds_relative_start + 1L
    cds_width[no_cds_idx] <- 0L

    C1_attribs <- recognizedBiomartAttribs(id_prefix)[["C1"]]
    has_C1_attribs <- all(C1_attribs %in% colnames(bm_result))
    C2_attribs <- recognizedBiomartAttribs(id_prefix)[["C2"]]
    has_C2_attribs <- all(C2_attribs %in% colnames(bm_result))

    if (has_C1_attribs)
        ans1 <- .extract_cds_ranges_from_C1(bm_result, id_prefix)
    if (has_C2_attribs) {
        ans2 <- .extract_cds_ranges_from_C2(bm_result, id_prefix)
        utr_anomaly <- mcols(ans2)$utr_anomaly
        tx_name_colname <- paste0(id_prefix, "transcript_id")
        tx_name <- bm_result[ , tx_name_colname]
        invalid_tx <- unique(tx_name[utr_anomaly != 0L])
        valid_tx_idx <- !(tx_name %in% invalid_tx)
        if (has_C1_attribs) {
            ## Check that 'ans1' agrees with 'ans2'.
            if (!identical(width(ans1)[valid_tx_idx],
                           width(ans2)[valid_tx_idx]))
                stop(wmsg("BioMart fatal data anomaly: ",
                          "CDS genomic coordinates are inconsistent with ",
                          "UTR genomic coordinates"))
            cds_idx <- which(valid_tx_idx & width(ans1) != 0L)
            if (!identical(start(ans1)[cds_idx], start(ans2)[cds_idx]))
                stop(wmsg("BioMart fatal data anomaly: ",
                          "CDS genomic coordinates are inconsistent with ",
                          "UTR genomic coordinates"))
        }
        ans1 <- ans2
    } else {
        valid_tx_idx <- seq_along(nrow(bm_result))
    }

    ## More checking of the CDS of the "valid" transcripts ("valid" here means
    ## with no UTR anomalies).
    .check_cds(ans1[valid_tx_idx], cds_width[valid_tx_idx],
               S4Vectors:::extract_data_frame_rows(bm_result, valid_tx_idx),
               id_prefix="ensembl_")
    ans1
}

.make_cds_df_from_ranges <- function(cds_ranges)
{
    no_cds_idx <- which(width(cds_ranges) == 0L)
    cds_start <- start(cds_ranges)
    cds_start[no_cds_idx] <- NA_integer_
    cds_end <- end(cds_ranges)
    cds_end[no_cds_idx] <- NA_integer_
    ans <- data.frame(cds_start=cds_start, cds_end=cds_end)
    utr_anomaly <- mcols(cds_ranges)$utr_anomaly
    if (!is.null(utr_anomaly))
        ans$utr_anomaly <- utr_anomaly
    ans
}

.makeBiomartSplicings <- function(filters, mart, transcripts_tx_id,
                                  recognized_attribs, id_prefix="ensembl_")
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
    available_attribs <- listAttributes(mart)$name
    has_group <- sapply(recognized_attribs[c("E2", "C1", "C2", "D1", "D2")],
                        function(attribs) all(attribs %in% available_attribs))
    get_groups <- c("E1", names(has_group)[has_group])
    attributes <- unlist(recognized_attribs[get_groups], use.names=FALSE)
    bm_result <- .getBM2(attributes, filters, mart=mart, bmHeader=FALSE)
    tx_name_colname <- paste0(id_prefix, "transcript_id")
    tx_name <- bm_result[[tx_name_colname]]
    splicings_tx_id <- transcripts_tx_id[tx_name]
    exon_name_colname <- paste0(id_prefix, "exon_id")
    splicings <- data.frame(
        tx_id=splicings_tx_id,
        exon_rank=bm_result$rank,
        exon_name=bm_result[[exon_name_colname]],
        exon_start=bm_result$exon_chrom_start,
        exon_end=bm_result$exon_chrom_end
    )
    if ((has_group[['C1']] || has_group[['C2']]) && has_group[['D1']]) {
        cds_ranges <- .extract_cds_ranges_from_bm_result(bm_result, id_prefix)
        splicings <- cbind(splicings, .make_cds_df_from_ranges(cds_ranges))
    }
    message("OK")
    splicings  
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'genes' data frame.
###

.makeBiomartGenes <- function(filters, mart,
                              transcripts_tx_id, recognized_attribs,
                              id_prefix="ensembl_")
{
    message("Download and preprocess the 'genes' data frame ... ",
            appendLF=FALSE)
    attributes <- c(recognized_attribs[['G']],
                    paste0(id_prefix, "transcript_id"))
    bm_result <- .getBM2(attributes, filters, mart=mart, bmHeader=FALSE)
    tx_name_colname <- paste0(id_prefix, "transcript_id")
    gene_id_colname <- paste0(id_prefix, "gene_id")
    tx_name <- bm_result[[tx_name_colname]]
    gene_id <- bm_result[[gene_id_colname]]
    keep_idx <- which(tx_name %in% names(transcripts_tx_id))
    tx_id <- transcripts_tx_id[tx_name[keep_idx]]
    gene_id <- gene_id[keep_idx]
    message("OK")
    data.frame(
        tx_id=tx_id,
        gene_id=gene_id
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Prepare the 'metadata' data frame.
###

.prepareBiomartMetadata <- function(mart, is_full_dataset, host, port,
                                    taxonomyId, miRBaseBuild)
{
    message("Prepare the 'metadata' data frame ... ",
            appendLF=FALSE)
    biomart <- biomaRt:::martBM(mart)
    dataset <- biomaRt:::martDataset(mart)
    mart_url <- biomaRt:::martHost(mart)
    mart_url <- sub("^[^/]+//", "", mart_url)
    mart_url <- unlist(strsplit(mart_url, "/"))[1]
    db_version <- .getBiomartDbVersion(mart, host, port, biomart)
    datasets <- listDatasets(mart)
    dataset_rowidx <- which(as.character(datasets$dataset) == dataset)
    ## This should never happen (the earlier call to useMart() would have
    ## failed in the first place).
    if (length(dataset_rowidx) != 1L)
        stop(wmsg("the BioMart database \"", biomaRt:::martBM(mart),
                  "\" has no (or more than one) \"", dataset, "\" datasets"))
    description <- as.character(datasets$description)[dataset_rowidx]
    dataset_version <- as.character(datasets$version)[dataset_rowidx]
    organism <- .extractOrganismFromDatasetDesc(description)
    if(is.na(taxonomyId)){
        taxonomyId <- GenomeInfoDb:::.taxonomyId(organism)
    }else{
        GenomeInfoDb:::.checkForAValidTaxonomyId(taxonomyId)
    }

    if (!isSingleStringOrNA(miRBaseBuild))
        stop(wmsg("'miRBaseBuild' must be a a single string or NA"))
    message("OK")
    data.frame(
        name=c("Data source",
               "Organism",
               "Taxonomy ID", 
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
                taxonomyId,
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
### makeTxDbFromBiomart()
###
### TODO: Rename the 'filters' arg -> 'filter' and set its default to NULL
### so it's consistent with the 'vals' arg of the transcripts(), exons(),
### cds(), and genes() extractors which we should also rename 'filter'.
###

makeTxDbFromBiomart <- function(biomart="ENSEMBL_MART_ENSEMBL",
                                dataset="hsapiens_gene_ensembl",
                                transcript_ids=NULL,
                                circ_seqs=DEFAULT_CIRC_SEQS,
                                filters="",
                                id_prefix="ensembl_",
                                host="www.ensembl.org",
                                port=80,
                                taxonomyId=NA,
                                miRBaseBuild=NA)
{
    mart <- .useMart2(biomart=biomart, dataset=dataset, host=host, port=port)
    id_prefix <- .normarg_id_prefix(id_prefix)
    filters <- .add_tx_id_filter(filters, transcript_ids, id_prefix)
    valid_filter_names <- listFilters(mart, what="name")
    invalid_filter_names <- setdiff(names(filters), valid_filter_names)
    if (length(invalid_filter_names) != 0L) {
        in1string <- paste0(invalid_filter_names, collapse=", ")
        stop(wmsg("Invalid filter name(s): ", in1string,
                  "\n\nPlease use the listFilters() function from the ",
                  "biomaRt package to get valid filter names."))
    }
    is_full_dataset <- length(filters) == 0L
    recognized_attribs <- recognizedBiomartAttribs(id_prefix)

    transcripts <- .makeBiomartTranscripts(filters, mart,
                                           transcript_ids,
                                           recognized_attribs,
                                           id_prefix)
    transcripts_tx_id <- transcripts$tx_id
    names(transcripts_tx_id) <- transcripts$tx_name
    chrominfo <- .makeBiomartChrominfo(mart,
                                       extra_seqnames=transcripts$tx_chrom,
                                       circ_seqs=circ_seqs,
                                       host, port)
    if (!is_full_dataset) {
        keep_idx <- which(chrominfo[ , "chrom"] %in% transcripts$tx_chrom)
        chrominfo <- S4Vectors:::extract_data_frame_rows(chrominfo, keep_idx)
    }
    splicings <- .makeBiomartSplicings(filters, mart,
                                       transcripts_tx_id,
                                       recognized_attribs,
                                       id_prefix=id_prefix)

    ## Drop transcripts with UTR anomalies.
    utr_anomaly <- splicings$utr_anomaly
    if (!is.null(utr_anomaly)) {
        invalid_tx <- unique(splicings[utr_anomaly != 0L, "tx_id"])
        if (length(invalid_tx) != 0L) {
            message("Drop transcripts with UTR anomalies (",
                    length(invalid_tx), " transcripts) ... ",
                    appendLF=FALSE)
            keep_idx1 <- !(transcripts$tx_id %in% invalid_tx)
            transcripts <- S4Vectors:::extract_data_frame_rows(transcripts,
                                                               keep_idx1)
            transcripts_tx_id <- transcripts_tx_id[keep_idx1]
            keep_idx2 <- !(splicings$tx_id %in% invalid_tx)
            splicings <- S4Vectors:::extract_data_frame_rows(splicings,
                                                             keep_idx2)
            message("OK")
        }
        splicings$utr_anomaly <- NULL
    }
        
    genes <- .makeBiomartGenes(filters, mart, transcripts_tx_id,
                               recognized_attribs, id_prefix)
    metadata <- .prepareBiomartMetadata(mart, is_full_dataset,
                                        host, port, taxonomyId, miRBaseBuild)

    message("Make the TxDb object ... ", appendLF=FALSE)
    txdb <- makeTxDb(transcripts, splicings,
                     genes=genes, chrominfo=chrominfo,
                     metadata=metadata, reassign.ids=TRUE)
    message("OK")
    txdb
}

