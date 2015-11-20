### =========================================================================
### Some non-exported tools to help explore/scan the BioMart landscape and
### find marts compatible with makeTxDbFromBiomart()
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Groups of BioMart attributes "recognized" by makeTxDbFromBiomart()
###
### REQUIRED attributes: makeTxDbFromBiomart() expects their presence
### in the dataset and will fail if they are not.
### OPTIONAL attributes: makeTxDbFromBiomart() will use them if
### they're present in the dataset but won't fail if they are not.
###

recognizedBiomartAttribs <- function(id_prefix="ensembl_")
{
    ans <- list(
        ## Group T: Transcript-level attributes
        ## ------------------------------------
        ## These attributes are REQUIRED.
        T=c(
            tx_name="<id_prefix>transcript_id",
            tx_type="transcript_biotype",
            tx_chrom="chromosome_name",
            tx_strand="strand",
            tx_start="transcript_start",
            tx_end="transcript_end"
        ),

        ## Groups E1 and E2: Exon-level attributes
        ## ---------------------------------------
        ## Attributes in group E1 are REQUIRED. Attributes in group E2 are
        ## OPTIONAL.
        E1=c(
            tx_name="<id_prefix>transcript_id",
            tx_name="strand",
            exon_rank="rank",
            exon_start="exon_chrom_start",
            exon_end="exon_chrom_end"
        ),
        E2=c(exon_name="<id_prefix>exon_id"),

        ## Groups C1 and C2: CDS-level attributes
        ## --------------------------------------
        ## These attributes are OPTIONAL.
        ##
        ## Ensembl added group C1 in release 74 (Dec 2013):
        ##   http://dec2013.archive.ensembl.org/info/website/news.html
        ## We added group C1 to recognizedBiomartAttribs() on 9/23/2014.
        ##
        ## Either group C1 or C2 is needed for inferring the CDS genomic
        ## coordinates (they cannot be inferred from group D1, see IMPORTANT
        ## NOTE ABOUT GROUP D1 below). If C1 and C2 are missing, the TxDb
        ## object will still be made but won't have any CDS (no row in the
        ## cds table).
        C1=c(
            cds_start="genomic_coding_start",
            cds_end="genomic_coding_end"
        ),
        ## Group C2 used to be group C. Was renamed C2 on 9/23/2014 when
        ## group C1 was added.
        C2=c(
            "5_utr_start",
            "5_utr_end",
            "3_utr_start",
            "3_utr_end"
        ),

        ## Groups D1 and D2: CDS-level attributes used for sanity checks only
        ## ------------------------------------------------------------------
        ## These attributes are OPTIONAL.
        ##
        ## IMPORTANT NOTE ABOUT GROUP D1: The "cds_start" and "cds_end"
        ## attributes that we get from BioMart are the CDS coordinates
        ## relative to the coding mRNA. This is *not* what we want to store
        ## in a TxDb object. What we want instead are the CDS *genomic*
        ## coordinates. Prior to Ensembl release 74, these were not provided
        ## by Ensembl so we were inferring them from the exon and UTR genomic
        ## coordinates (the UTR genomic coordinates are obtained thru the
        ## attributes in group C2).
        ## So attributes in group D1 (as well as in group D2) are only used
        ## to perform sanity checks.
        D1=c(
            "cds_start",
            "cds_end"
        ),
        D2="cds_length",

        ## Group G: Gene-level attributes
        ## ------------------------------
        ## These attributes are REQUIRED.
        G=c(gene_id="<id_prefix>gene_id")
    )
    lapply(ans, gsub, pattern="<id_prefix>", replacement=id_prefix, fixed=TRUE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Find marts compatible with makeTxDbFromBiomart().
###

### 'mart' can be either a Mart object or the name of a Mart service (single
### string). Returns a named list of 2-col data frames with one elt per
### dataset in 'mart'. Each data frame describes the attributes that are
### available for the corresponding dataset.
### Typical use:
###   ensembl_attribs_list <- GenomicFeatures:::.get_attribs_for_each_dataset(
###                               "ENSEMBL_MART_ENSEMBL")
###   sapply(ensembl_attribs_list, nrow)
###   table(sapply(ensembl_attribs_list, nrow))
### or to walk only on the first 5 datasets:
###   GenomicFeatures:::.get_attribs_for_each_dataset("ENSEMBL_MART_ENSEMBL",
###                               length.out=5)
.get_attribs_for_each_dataset <- function(mart, length.out=NULL)
{
    if (!is(mart, "Mart"))
        mart <- useMart(mart)
    available_datasets <- listDatasets(mart)
    if (!is.null(length.out))
        available_datasets <- head(available_datasets, n=length.out)
    ans_len <- nrow(available_datasets)
    ans <- vector(mode="list", length=ans_len)
    names(ans) <- as.character(available_datasets$dataset)
    for (i in seq_len(ans_len)) {
        dataset <- names(ans)[i]
        mart <- useDataset(dataset, mart=mart)
        message("Getting attributes for dataset \"", dataset, "\" ... ",
                appendLF=FALSE)
        ans[[i]] <- listAttributes(mart)
        message("OK")
    }
    ans
}

### 'available_attribs' can be a character vector of attribute names, a 2-col
### data frame as returned by biomaRt::listAttributes(), or a Mart object
### with a selected dataset.
.find_available_attrib_groups <- function(available_attribs)
{
    if (is(available_attribs, "Mart"))
        available_attribs <- listAttributes(available_attribs)
    if (is.data.frame(available_attribs)) {
        stopifnot(identical(colnames(available_attribs),
                            c("name", "description")))
        available_attribs <- available_attribs$name
    }
    if (!is.character(available_attribs))
        stop("invalid 'available_attribs' argument")
    recognized_attribs <- recognizedBiomartAttribs()
    has_group <- sapply(recognized_attribs,
                        function(attribs) all(attribs %in% available_attribs))
    ans <- paste0(names(has_group)[has_group], collapse="+") 
    if (ans == "")
        ans <- "none"
    ans
}

### Vectorized version of .find_available_attrib_groups().
### 'available_attribs' can be a named list (as returned by
### .get_attribs_for_each_dataset()), a Mart object, or the name of a Mart
### service (single string).
### Typical use:
###   ensembl_attrib_groups <-
###     GenomicFeatures:::.find_available_attrib_groups_for_each_dataset(
###                           "ENSEMBL_MART_ENSEMBL")
.find_available_attrib_groups_for_each_dataset <-
    function(available_attribs, length.out=NULL)
{
    if (!is.list(available_attribs))
        available_attribs <- .get_attribs_for_each_dataset(
                                     available_attribs,
                                     length.out=length.out)
    sapply(available_attribs, .find_available_attrib_groups)
}

### 'biomart' and 'version' must be single character strings.
### Typical use:
###   library(biomaRt)
###   marts <- listMarts()
###   GenomicFeatures:::scanMart(marts$biomart[1], marts$version[1],
###                              length.out=5)
scanMart <- function(biomart, version, length.out=NULL)
{
    biomart <- as.character(biomart)
    version <- as.character(version)
    cat("Scanning ", biomart, " mart ... ", sep="")
    suppressMessages(
      available_groups <- .find_available_attrib_groups_for_each_dataset(
                               biomart,
                               length.out=length.out)
    )
    cat("OK\n")
    cat("-------------------------- scan result --------------------------\n")
    cat("### biomart: ", biomart, "\n", sep="")
    cat("### version: ", version, "\n", sep="")
    available_datasets <- names(available_groups)
    if (length(available_datasets) > 4L)
        available_datasets <- c(available_datasets[1:3], "...")
    cat("### nb of available datasets: ", length(available_groups),
        " (", paste(available_datasets, collapse=", "), ")\n",
        sep="")
    cat("### summary of attribute groups:\n")
    tbl <- table(available_groups)
    headers <- c("available attribute groups", "nb of datasets")
    col1 <- c(headers[1L], format(names(tbl), width=nchar(headers[1L])))
    col2 <- c(headers[2L], format(as.integer(tbl), width=nchar(headers[2L])))
    cat(paste0("###     ", col1, " | ", col2, "\n"), sep="")
    cat("\n")
}

findCompatibleMarts <- function(marts=NULL, ...)
{
    if (is.null(marts))
        marts <- listMarts(...)
    biomarts <- as.character(marts$biomart)
    versions <- as.character(marts$version)
    for (i in seq_len(nrow(marts))) {
        res <- try(scanMart(biomarts[i], versions[i]), silent=TRUE)
        if (is(res, "try-error"))
            message("scanMart() failed! ==> skipping ", biomarts[i], " mart\n")
    }
}


### =========================================================================
### findCompatibleMarts() output as of 9/24/2014
###
### Using GenomicFeatures 1.17.15 (rev 94513).
###
### Notes:
###  - Output was manually cleaned to keep only biomarts providing at least
###    1 dataset with the REQUIRED attributes (i.e. groups T, E1, and G).
###  - Since Ensembl Genomes release 17 (Feb 2013), BioMart access to
###    Ensembl Bacteria is no longer possible. See announcement:
###      http://www.ensembl.info/blog/2013/02/07/ensembl-genomes-release-17/
### -------------------------------------------------------------------------
###
### biomart: ensembl
### version: ENSEMBL GENES 75 (SANGER UK)
### nb of available datasets: 66 (hsapiens_gene_ensembl,
###     oanatinus_gene_ensembl, cporcellus_gene_ensembl,
###     gaculeatus_gene_ensembl, ...)
### summary of attribute groups:
###     available attribute groups | nb of datasets
###     T+E1+E2+C1+C2+D1+D2+G      |             66
###
### biomart: fungi_mart_22
### version: ENSEMBL FUNGI 22 (EBI UK)
### nb of available datasets: 45 (aterreus_eg_gene,
###     mlaricipopulina_eg_gene, treesei_eg_gene, ...)
### summary of attribute groups:
###     available attribute groups | nb of datasets
###     T+E1+E2+C1+C2+D1+G         |             45
###
### biomart: metazoa_mart_22
### version: ENSEMBL METAZOA 22 (EBI UK)
### nb of available datasets: 52 (tcastaneum_eg_gene, dgrimshawi_eg_gene,
###     rprolixus_eg_gene, ...)
### summary of attribute groups:
###     available attribute groups | nb of datasets
###     T+E1+E2+C1+C2+D1+G         |             52
###
### biomart: plants_mart_22
### version: ENSEMBL PLANTS 22 (EBI UK)
### nb of available datasets: 33 (athaliana_eg_gene, atauschii_eg_gene,
###     obrachyantha_eg_gene, ppersica_eg_gene, ...)
### summary of attribute groups:
###     available attribute groups | nb of datasets
###     T+E1+E2+C1+C2+D1+G         |             33
###
### biomart: protists_mart_22
### version: ENSEMBL PROTISTS 22 (EBI UK)
### nb of available datasets: 29 (pramorum_eg_gene,
###     pvivax_eg_gene, piwayamai_eg_gene,  ... )
### summary of attribute groups:
###     available attribute groups | nb of datasets
###     T+E1+E2+C1+C2+D1+G         |             29
###
### biomart: Breast_mart_69
### version: BCCTB Bioinformatics Portal (UK and Ireland)
### nb of available datasets: 2 (breastCancer_expressionStudy,
###     hsapiens_gene_breastCancer)
### summary of attribute groups:
###     available attribute groups | nb of datasets
###     none                       |              1
###     T+E1+E2+C1+C2+D1+D2+G      |              1
###
### biomart: Pancreas63
### version: PANCREATIC EXPRESSION DATABASE (BARTS CANCER INSTITUTE UK)
### nb of available datasets: 2 (hsapiens_gene_pancreas,
###     hsapiens_cancer_pancreas)
### summary of attribute groups:
###     available attribute groups | nb of datasets
###     none                       |              1
###     T+E1+E2+C1+C2+D1+D2+G      |              1
###
### biomart: ENSEMBL_MART_PLANT
### version: GRAMENE 40 ENSEMBL GENES (CSHL/CORNELL US)
### nb of available datasets: 33 (atauschii_eg_gene, obrachyantha_eg_gene,
###     ptrichocarpa_eg_gene,  ... )
### summary of attribute groups:
###     available attribute groups | nb of datasets
###     T+E1+E2+C1+C2+D1+G         |             33


### =========================================================================
### findCompatibleMarts() output as of 6/28/2010
###
### Notes:
###  - Group A is the union of groups T and E1.
###  - Only biomarts containing datasets with at least attribute groups A and
###    G are listed.
### -------------------------------------------------------------------------
###
### biomart: ensembl
### version: ENSEMBL GENES 58 (SANGER UK)
### nb of datasets: 51 (hsapiens_gene_ensembl, oanatinus_gene_ensembl,
###                     tguttata_gene_ensembl, cporcellus_gene_ensembl, ...)
### NOTE: the mgallopavo_gene_ensembl dataset seems to be broken!
### table of attribute groups: AE2CDG:50
###
### biomart: bacterial_mart_5
### version: ENSEMBL BACTERIA 5 (EBI UK)
### nb of datasets: 183 (str_57_gene, esc_20_gene, myc_25994_gene, ...)
### table of attribute groups: AE2G:183
###
### biomart: fungal_mart_5
### version: ENSEMBL FUNGAL 5 (EBI UK)
### nb of datasets: 12 (aniger_eg_gene, aflavus_eg_gene, aterreus_eg_gene, ...)
### table of attribute groups: AE2G:12 
###
### biomart: metazoa_mart_5
### version: ENSEMBL METAZOA 5 (EBI UK)
### nb of datasets: 23 (dgrimshawi_eg_gene, ppacificus_eg_gene,
###                     dpseudoobscura_eg_gene, ...)
### table of attribute groups: AE2G:23
###
### biomart: plant_mart_5
### version: ENSEMBL PLANT 5 (EBI UK)
### nb of datasets: 8 (sbicolor_eg_gene, bdistachyon_eg_gene,
###                    alyrata_eg_gene, ...)
### table of attribute groups: AE2G:8 
###
### biomart: protist_mart_5
### version: ENSEMBL PROTISTS 5 (EBI UK)
### nb of datasets: 6 (tpseudonana_gene, ptricornutum_gene, pknowlesi_gene, ...)
### table of attribute groups: AE2G:6
###
### biomart: ensembl_expressionmart_48
### version: EURATMART (EBI UK)
### nb of datasets: 1 (rnorvegicus_expr_gene_ensembl)
### table of attribute groups: AG:1
###
### biomart: Ensembl56
### version: PANCREATIC EXPRESSION DATABASE (INSTITUTE OF CANCER UK)
### nb of datasets: 1 (hsapiens_gene_pancreas)
### table of attribute groups: AE2CDG:1
###
### biomart: ENSEMBL_MART_ENSEMBL
### version: GRAMENE 30 ENSEMBL GENES (CSHL/CORNELL US)
### nb of datasets: 8 (sbicolor_eg_gene, bdistachyon_eg_gene,
###                    alyrata_eg_gene, ...)
### table of attribute groups: AE2G:8

