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


.getBiomartDbVersion <- function(biomart)
{
    marts <- listMarts()
    mart_rowidx <- which(as.character(marts$biomart) == biomart)
    ## This should never happen.
    if (length(mart_rowidx) != 1L)
        stop("found 0 or more than 1 \"", biomart, "\" BioMart database")
    as.character(marts$version)[mart_rowidx]
}

.extractEnsemblReleaseFromDbVersion <- function(db_version)
    sub("^ENSEMBL GENES ([^[:space:]]+) \\(SANGER UK\\)", "\\1", db_version)

### Groups of BioMart attributes:
###   - A1, A2 and G are required attributes;
###   - B, C and D are optional attributes: C is required for inferring the
###     CDS (they cannot be inferred from D). Therefore, if C is missing,
###     the TranscriptDb object can still be made but won't have any CDS (no
###     row in the cds table). D is only used for sanity check.
.A1_ATTRIBS <- c("ensembl_transcript_id",
                 "chromosome_name",
                 "strand",
                 "transcript_start",
                 "transcript_end")

.A2_ATTRIBS <- c("ensembl_transcript_id",
                 "strand",
                 "rank",
                 "exon_chrom_start",
                 "exon_chrom_end")

.B_ATTRIB <- "ensembl_exon_id"

.C_ATTRIBS <- c("5_utr_start",
                "5_utr_end",
                "3_utr_start",
                "3_utr_end")

.D_ATTRIBS <- c("cds_start",
                "cds_end",
                "cds_length")

.G_ATTRIB <- "ensembl_gene_id"

### 'attribs' can be either a Mart object or a 2-col data frame as returned by
### 'listAttributes()'.
.getDatasetAttrGroups <- function(attribs)
{
    if (is(attribs, "Mart"))
        attribs <- listAttributes(attribs)
    else if (!is.data.frame(attribs) ||
             !identical(colnames(attribs), c("name", "description")))
        stop("invalid 'attribs' object")
    attrgroups <- "none"
    ## Group A: Required attributes.
    attrA <- unique(c(.A1_ATTRIBS, .A2_ATTRIBS))
    if (all(attrA %in% attribs$name))
        attrgroups <- "A"
    ## Groups B, C and D are optional attributes.
    ## C is required for inferring the CDS (they cannot be inferred from D).
    ## Therefore, if C is missing, the TranscriptDb object can still be made
    ## but won't have any CDS (no row in the cds table).
    if (.B_ATTRIB %in% attribs$name)
        attrgroups <- paste(attrgroups, "B", sep="")
    if (all(.C_ATTRIBS %in% attribs$name))
        attrgroups <- paste(attrgroups, "C", sep="")
    if (all(.D_ATTRIBS %in% attribs$name))
        attrgroups <- paste(attrgroups, "D", sep="")
    ## Group G: Required attribute.
    if (.G_ATTRIB %in% attribs$name)
        attrgroups <- paste(attrgroups, "G", sep="")
    attrgroups
}

### 'attrlist' can be a list (as returned by getMartAttribList()), a Mart
### object, or the name of a Mart service (single string).
### Typical use:
###   ensembl_attrgroups <-
###       GenomicFeatures:::.getAllDatasetAttrGroups("ensembl")
.getAllDatasetAttrGroups <- function(attrlist)
{
    if (!is.list(attrlist))
        attrlist <- getMartAttribList(attrlist)
    sapply(attrlist, .getDatasetAttrGroups)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'transcripts' data frame.
###

.makeBiomartTranscripts <- function(filters, values, mart, transcript_ids)
{
    message("Download and preprocess the 'transcripts' data frame ... ",
            appendLF=FALSE)
    bm_table <- getBM(.A1_ATTRIBS, filters=filters, values=values, mart=mart)
    if (!is.null(transcript_ids)) {
        idx <- !(transcript_ids %in% bm_table$ensembl_transcript_id)
        if (any(idx)) {
            bad_ids <- transcript_ids[idx]
            stop("invalid transcript ids: ",
                 paste(bad_ids, collapse=", "), sep="")
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
    transcripts_tx_id <- seq_len(nrow(bm_table))
    transcripts_tx_name <- bm_table$ensembl_transcript_id
    if (any(duplicated(transcripts_tx_name)))
        stop("the 'ensembl_transcript_id' attribute contains duplicated values")
    transcripts <- data.frame(
        tx_id=transcripts_tx_id,
        tx_name=transcripts_tx_name,
        tx_chrom=as.character(bm_table$chromosome_name),
        tx_strand=ifelse(bm_table$strand == 1, "+", "-"),
        tx_start=bm_table$transcript_start,
        tx_end=bm_table$transcript_end
    )
    message("OK")
    transcripts
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'chrominfo' data frame.
###

### Uses RCurl to access and list the content of an FTP dir.
.lsFtpUrl <- function(url)
{
    doc <- getURL(url)
    listing <- strsplit(doc, "\n", fixed=TRUE)[[1L]]
    ## Keep field no. 8 only
    pattern <- paste(c("^", rep.int("[^[:space:]]+[[:space:]]+", 8L)),
                     collapse="")
    listing <- sub(pattern, "", listing)
    sub("[[:space:]].*$", "", listing)
}

.ENSEMBL.PUB_FTP_URL <- "ftp://ftp.ensembl.org/pub/"

.Ensembl.getFtpUrlToMySQL <- function(release=NA)
{
    if (is.na(release))
        pub_subdir <- "current_mysql"
    else
        pub_subdir <- paste("release-", release, "/mysql", sep="")
    paste(.ENSEMBL.PUB_FTP_URL, pub_subdir, "/", sep="")
}

.Ensembl.listMySQLCoreDirs <- function(release=NA, url=NA)
{
    if (is.na(url))
        url <- .Ensembl.getFtpUrlToMySQL(release)
    core_dirs <- .lsFtpUrl(url)
    pattern <- "_core_"
    if (!is.na(release))
        pattern <- paste(pattern, release, "_", sep="")
    core_dirs[grep(pattern, core_dirs, fixed=TRUE)]
}

.Ensembl.getMySQLCoreDir <- function(dataset, release=NA, url=NA)
{
    if (is.na(url))
        url <- .Ensembl.getFtpUrlToMySQL(release)
    core_dirs <- .Ensembl.listMySQLCoreDirs(release=release, url=url)
    shortnames <- sapply(strsplit(core_dirs, "_", fixed=TRUE),
                         function(x)
                           paste(substr(x[1L], 1L, 1L), x[2L], sep=""))
    shortname0 <- strsplit(dataset, "_", fixed=TRUE)[[1L]][1L]
    core_dir <- core_dirs[shortnames == shortname0]
    if (length(core_dir) != 1L)
        stop("found 0 or more than 1 subdir for \"", dataset,
             "\" dataset at ", url)
    core_dir
}

.Ensembl.getMySQLCoreUrl <- function(dataset, release=NA, url=NA)
{
    if (is.na(url))
        url <- .Ensembl.getFtpUrlToMySQL(release)
    core_dir <- .Ensembl.getMySQLCoreDir(dataset, release=release, url=url)
    paste(url, core_dir, "/", sep="")
}

.Ensembl.fetchTableDump <- function(base_url, tablename, colnames)
{
    url <- paste(base_url, tablename, ".txt.gz", sep="")
    destfile <- tempfile()
    download.file(url, destfile, quiet=TRUE)
    data <- read.table(destfile, sep="\t", quote="",
                       col.names=colnames, comment.char="",
                       stringsAsFactors=FALSE)
    unlink(destfile)
    data
}

.Ensembl.fetchAttribTypeIdForTopLevelSequence <- function(core_url)
{
    ## Get 'attrib_type' table.
    colnames <- c("attrib_type_id", "code", "name", "description")
    attrib_type <- .Ensembl.fetchTableDump(core_url, "attrib_type", colnames)
    i <- which(attrib_type$code == "toplevel")
    if (length(i) != 1L)
        stop("Ensembl data anomaly: \"toplevel\" attrib found 0 or more ",
             "than once in attrib_type table at ", core_url)
    attrib_type$attrib_type_id[i]
}

.Ensembl.fetchTopLevelSequenceIds <- function(core_url)
{
    id0 <- .Ensembl.fetchAttribTypeIdForTopLevelSequence(core_url)
    ## Get 'seq_region_attrib' table.
    colnames <- c("seq_region_id", "attrib_type_id", "value")
    seq_region_attrib <- .Ensembl.fetchTableDump(core_url,
                             "seq_region_attrib", colnames)
    seq_region_attrib$seq_region_id[seq_region_attrib$attrib_type_id == id0]
}

### Fetch names and lengths for the "toplevel" sequences in the
### 'seq_region' table of the Ensemble Core DB specified by 'core_url'.
### Ensembl Core Schema Documentation:
###   http://www.ensembl.org/info/docs/api/core/core_schema.html
### The full schema:
###   ftp://ftp.ensembl.org/pub/ensembl/sql/table.sql
### Typical use:
###   core_url <- .Ensembl.getMySQLCoreUrl("hsapiens_gene_ensembl")
###   extra_seqnames <- unique(as.character(transcripts$tx_chrom))
###   extra_seqnames <- c("GL000217.1", "NC_012920")
###   .fetchChromLengthsFromCoreUrl(core_url, extra_seqnames=extra_seqnames)
.fetchChromLengthsFromCoreUrl <- function(core_url, extra_seqnames=NULL)
{
    ## Get 'seq_region' table.
    colnames <- c("seq_region_id", "name", "coord_system_id", "length")
    seq_region <- .Ensembl.fetchTableDump(core_url, "seq_region", colnames)
    ## Get 'coord_system' table.
    colnames <- c("coord_system_id", "species_id", "name",
                  "version", "rank", "attrib")
    coord_system <- .Ensembl.fetchTableDump(core_url, "coord_system", colnames)
    ## 1st filtering: keep only "default_version" sequences.
    idx1 <- grep("default_version", coord_system$attrib, fixed=TRUE)
    ids1 <- coord_system$coord_system_id[idx1]
    seq_region <- seq_region[seq_region$coord_system_id %in% ids1, , drop=FALSE]
    ## 2nd filtering: keep "toplevel" sequences + extra sequences.
    seq_region_ids <- .Ensembl.fetchTopLevelSequenceIds(core_url)
    ans_idx <- seq_region$seq_region_id %in% seq_region_ids
    if (!is.null(extra_seqnames)) {
        if (!all(extra_seqnames %in% seq_region$name))
            stop("failed to fetch all chromosome lengths")
        ## Add extra sequences to the index.
        ans_idx <- ans_idx | (seq_region$name %in% extra_seqnames)
    }
    ans <- seq_region[ans_idx, c("name", "length"), drop=FALSE]
    row.names(ans) <- NULL
    if (any(duplicated(ans$name)))
        stop("failed to fetch all chromosome lengths unambiguously")
    ans
}

.fetchChromLengthsFromEnsembl <- function(dataset, release=NA,
                                          extra_seqnames=NULL)
{
    core_url <- .Ensembl.getMySQLCoreUrl(dataset, release=release)
    .fetchChromLengthsFromCoreUrl(core_url, extra_seqnames=extra_seqnames)
}

### Returns NULL if it fails to fetch the chromosome lengths from the
### remote resource.
.makeBiomartChrominfo <- function(mart, circ_seqs, extra_seqnames=NULL)
{
    biomart <- biomaRt:::martBM(mart)
    dataset <- biomaRt:::martDataset(mart)
    if (biomart == "ensembl") {
        message("Download and preprocess the 'chrominfo' data frame ... ",
                appendLF=FALSE)
        db_version <- .getBiomartDbVersion(biomart)
        ensembl_release <- .extractEnsemblReleaseFromDbVersion(db_version)
        chromlengths <- try(.fetchChromLengthsFromEnsembl(dataset,
                                release=ensembl_release,
                                extra_seqnames=extra_seqnames),
                            silent=TRUE)
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

discoverBiomartChrominfo <- function(biomart="ensembl",
                                     dataset="hsapiens_gene_ensembl")
{
    mart <- .parseBMMartParams(biomart=biomart,
                              dataset=dataset)
    filters <- .parseBMFiltersParams(transcript_ids=NULL)
    values <- .parseBMValuesParams(transcript_ids=NULL)        
    transcripts <- .makeBiomartTranscripts(filters, values, mart,
                                           transcript_ids=NULL)
    chrominfo <- .makeBiomartChrominfo(mart, circ_seqs=character(),
                                       extra_seqnames=transcripts$tx_chrom)
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

.extractCdsRangesFromBiomartTable <- function(bm_table)
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

.makeCdsDataFrameFromRanges <- function(cds_ranges)
{
    nocds_idx <- width(cds_ranges) == 0L
    cds_start <- start(cds_ranges)
    cds_start[nocds_idx] <- NA_integer_
    cds_end <- end(cds_ranges)
    cds_end[nocds_idx] <- NA_integer_
    data.frame(cds_start=cds_start, cds_end=cds_end)
}

### Ironically the cds_start and cds_end attributes that we get from
### BioMart are pretty useless because they are relative to the coding
### mRNA. However, the utr coordinates are relative to the chromosome so
### we use them to infer the cds coordinates. We also retrieve the
### cds_length attribute as a sanity check.
.makeBiomartSplicings <- function(filters, values, mart, transcripts_tx_name)
{
    ## Those are the strictly required fields.
    splicings0 <- data.frame(
        tx_id=integer(0),
        exon_rank=integer(0),
        exon_start=integer(0),
        exon_end=integer(0)
    )
    if (length(transcripts_tx_name) == 0L)
        return(splicings0)
    message("Download and preprocess the 'splicings' data frame ... ",
            appendLF=FALSE)
    allattribs <- listAttributes(mart)$name
    attributes <- .A2_ATTRIBS
    if (.B_ATTRIB %in% allattribs)
        attributes <- c(attributes, .B_ATTRIB)
    if (all(.C_ATTRIBS %in% allattribs))
        attributes <- c(attributes, .C_ATTRIBS)
    if ("cds_length" %in% allattribs)
        attributes <- c(attributes, "cds_length")
    bm_table <- getBM(attributes, filters=filters, values=values, mart=mart)
    splicings_tx_id <- as.integer(factor(bm_table$ensembl_transcript_id,
                                         levels=transcripts_tx_name))
    splicings <- data.frame(
        tx_id=splicings_tx_id,
        exon_rank=bm_table$rank,
        exon_name=bm_table$ensembl_exon_id,
        exon_start=bm_table$exon_chrom_start,
        exon_end=bm_table$exon_chrom_end
    )
    if (all(.C_ATTRIBS %in% allattribs) && ("cds_length" %in% allattribs)) {
        cds_ranges <- .extractCdsRangesFromBiomartTable(bm_table)
        splicings <- cbind(splicings, .makeCdsDataFrameFromRanges(cds_ranges))
    }
    message("OK")
    splicings
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'genes' data frame.
###

.makeBiomartGenes <- function(filters, values, mart, transcripts_tx_name)
{
    message("Download and preprocess the 'genes' data frame ... ",
            appendLF=FALSE)
    attributes <- c(.G_ATTRIB, "ensembl_transcript_id")
    bm_table <- getBM(attributes, filters=filters, values=values, mart=mart)
    genes_tx_id <- as.integer(factor(bm_table$ensembl_transcript_id,
                                     levels=transcripts_tx_name))
    message("OK")
    data.frame(
        tx_id=genes_tx_id,
        gene_id=bm_table$ensembl_gene_id
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Prepare the 'metadata' data frame.
###

.prepareBiomartMetadata <- function(mart, is_full_dataset)
{
    message("Prepare the 'metadata' data frame ... ",
            appendLF=FALSE)
    biomart <- biomaRt:::martBM(mart)
    dataset <- biomaRt:::martDataset(mart)
    db_version <- .getBiomartDbVersion(biomart)
    datasets <- listDatasets(mart)
    dataset_rowidx <- which(as.character(datasets$dataset) == dataset)
    ## This should never happen (the above call to useMart() would have failed
    ## in the first place).
    if (length(dataset_rowidx) != 1L)
        stop("the BioMart database \"", biomaRt:::martBM(mart),
             "\" has no (or more than one) \"", dataset, "\" datasets")
    description <- as.character(datasets$description)[dataset_rowidx]
    dataset_version <- as.character(datasets$version)[dataset_rowidx]
    message("OK")
    data.frame(
        name=c("Data source",
               "BioMart database",
               "BioMart database version",
               "BioMart dataset",
               "BioMart dataset description",
               "BioMart dataset version",
               "Full dataset"),
        value=c("BioMart",
                biomart,
                db_version,
                dataset,
                description,
                dataset_version,
                ifelse(is_full_dataset, "yes", "no"))
    )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTranscriptDbFromBiomart()
###

.parseBMMartParams <- function(biomart="ensembl",
                                      dataset="hsapiens_gene_ensembl")
{
    if (is.factor(biomart))
        biomart <- as.character(biomart)
    if (is(dataset, "AsIs"))
        dataset <- as.character(dataset)
    if (!isSingleString(biomart))
        stop("'biomart' must be a single string")
    useMart(biomart=biomart, dataset=dataset)
}

.parseBMFiltersParams <- function(transcript_ids)
{
    if (is.null(transcript_ids)) {
        filters <- ""
    } else if (is.character(transcript_ids)
            && !any(is.na(transcript_ids))) {
        filters <- "ensembl_transcript_id"
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
##                                       circ_seqs=DEFAULTCIRCSTRS,
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
                                        circ_seqs=DEFAULTCIRCSTRS,
                                        transcript_ids=NULL)
{
    ## Could be that the user got the 'biomart' and/or 'dataset' values
    ## programmatically via calls to listMarts() and/or listDatasets().
    mart <- .parseBMMartParams(biomart=biomart,
                              dataset=dataset)
    filters <- .parseBMFiltersParams(transcript_ids)
    values <- .parseBMValuesParams(transcript_ids)
    
    transcripts <- .makeBiomartTranscripts(filters, values, mart,
                                           transcript_ids)
    chrominfo <- .makeBiomartChrominfo(mart, circ_seqs=circ_seqs,
                                       extra_seqnames=transcripts$tx_chrom)
    splicings <- .makeBiomartSplicings(filters, values, mart,
                                       transcripts$tx_name)
    genes <- .makeBiomartGenes(filters, values, mart, transcripts$tx_name)
    metadata <- .prepareBiomartMetadata(mart, is.null(transcript_ids))

    message("Make the TranscriptDb object ... ", appendLF=FALSE)
    txdb <- makeTranscriptDb(transcripts, splicings,
                             genes=genes, chrominfo=chrominfo,
                             metadata=metadata)
    message("OK")
    txdb
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some non-exported tools to help exploring/scanning the BioMart landscape.
###

### 'mart' can be either a Mart object or the name of a Mart service (single
### string). Returns a named list of 2-col data frames with one elt per
### dataset in 'mart'. Each data frame describes the attributes that are
### available for the corresponding dataset.
### Typical use:
###   ensembl_attrlist <- GenomicFeatures:::getMartAttribList("ensembl")
###   sapply(ensembl_attrlist, nrow)
getMartAttribList <- function(mart)
{
    if (!is(mart, "Mart"))
        mart <- useMart(mart)
    datasets <- listDatasets(mart)
    ans_length <- nrow(datasets)
    ans <- vector(mode="list", length=ans_length)
    names(ans) <- as.character(datasets$dataset)
    for (i in seq_len(ans_length)) {
        dataset <- names(ans)[i]
        mart <- useDataset(dataset, mart=mart)
        message("Getting attributes for dataset \"", dataset, "\"... ",
                appendLF=FALSE)
        ans[[i]] <- listAttributes(mart)
        message("OK")
    }
    ans
}

### 'biomart' and 'version' must be single character strings.
scanMart <- function(biomart, version)
{
    cat("Scanning ", biomart, "... ", sep="")
    suppressMessages(attrgroups <- .getAllDatasetAttrGroups(biomart))
    cat("OK\n")
    cat("biomart: ", biomart, "\n", sep="")
    cat("version: ", version, "\n", sep="")
    tmp <- names(attrgroups)
    if (length(tmp) > 3L)
        tmp <- c(tmp[1:3], "...")
    cat("nb of datasets: ", length(attrgroups),
        " (", paste(tmp, collapse=", "), ")\n",
        sep="")
    if (length(attrgroups) != 0L) {
        tbl <- table(attrgroups)
        tbl2 <- as.integer(tbl)
        names(tbl2) <- names(tbl)
        tmp <- paste(names(tbl2), ":", tbl2, sep="", collapse=", ")
        cat("table of attribute groups: ", tmp, "\n", sep="")
    }
    cat("\n")
}

scanMarts <- function(marts=NULL)
{
    if (is.null(marts))
        marts <- listMarts()
    biomarts <- as.character(marts$biomart)
    versions <- as.character(marts$version)
    for (i in seq_len(nrow(marts)))
        scanMart(biomarts[i], versions[i])
}

### scanMarts() output as of 6/28/2010 (only biomarts with at least groups
### A and G are listed):
###
### biomart: ensembl
### version: ENSEMBL GENES 58 (SANGER UK)
### nb of datasets: 51 (hsapiens_gene_ensembl, oanatinus_gene_ensembl,
###                     tguttata_gene_ensembl, cporcellus_gene_ensembl, ...)
### NOTE: the mgallopavo_gene_ensembl dataset seems to be broken!
### table of attribute groups: ABCDG:50
###
### biomart: bacterial_mart_5
### version: ENSEMBL BACTERIA 5 (EBI UK)
### nb of datasets: 183 (str_57_gene, esc_20_gene, myc_25994_gene, ...)
### table of attribute groups: ABG:183
###
### biomart: fungal_mart_5
### version: ENSEMBL FUNGAL 5 (EBI UK)
### nb of datasets: 12 (aniger_eg_gene, aflavus_eg_gene, aterreus_eg_gene, ...)
### table of attribute groups: ABG:12 
###
### biomart: metazoa_mart_5
### version: ENSEMBL METAZOA 5 (EBI UK)
### nb of datasets: 23 (dgrimshawi_eg_gene, ppacificus_eg_gene,
###                     dpseudoobscura_eg_gene, ...)
### table of attribute groups: ABG:23
###
### biomart: plant_mart_5
### version: ENSEMBL PLANT 5 (EBI UK)
### nb of datasets: 8 (sbicolor_eg_gene, bdistachyon_eg_gene,
###                    alyrata_eg_gene, ...)
### table of attribute groups: ABG:8 
###
### biomart: protist_mart_5
### version: ENSEMBL PROTISTS 5 (EBI UK)
### nb of datasets: 6 (tpseudonana_gene, ptricornutum_gene, pknowlesi_gene, ...)
### table of attribute groups: ABG:6
###
### biomart: ensembl_expressionmart_48
### version: EURATMART (EBI UK)
### nb of datasets: 1 (rnorvegicus_expr_gene_ensembl)
### table of attribute groups: AG:1
###
### biomart: Ensembl56
### version: PANCREATIC EXPRESSION DATABASE (INSTITUTE OF CANCER UK)
### nb of datasets: 1 (hsapiens_gene_pancreas)
### table of attribute groups: ABCDG:1
###
### biomart: ENSEMBL_MART_ENSEMBL
### version: GRAMENE 30 ENSEMBL GENES (CSHL/CORNELL US)
### nb of datasets: 8 (sbicolor_eg_gene, bdistachyon_eg_gene,
###                    alyrata_eg_gene, ...)
### table of attribute groups: ABG:8

