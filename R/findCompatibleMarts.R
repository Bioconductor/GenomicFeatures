### =========================================================================
### Some non-exported tools to help exploring/scanning the BioMart landscape
### and finding marts compatible with makeTranscriptDbFromBiomart()
### -------------------------------------------------------------------------
###

### Groups of BioMart attributes used by makeTranscriptDbFromBiomart():
###
###   - Group A (i.e. A1 + A2) are required attributes.
###
###   - Groups B, C1, C2, and D are optional attributes.
###     Either group C1 or C2 is required for inferring the CDS (they cannot
###     be inferred from D). If C1 and C2 are missing, the TxDb object will
###     still be made but won't have any CDS (no row in the cds table).
###     Attributes in group D are only used for sanity checks.
###
###   - Group G are required attributes.
###

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

### Added on 9/23/2014.
.C1_ATTRIBS <- c("genomic_coding_start",
                 "genomic_coding_end")

### C2 used to be group C. Was renamed C2 on 9/23/2014.
.C2_ATTRIBS <- c("5_utr_start",
                 "5_utr_end",
                 "3_utr_start",
                 "3_utr_end")

.D_ATTRIBS <- c("cds_start",
                "cds_end",
                "cds_length")

.G_ATTRIB <- "ensembl_gene_id"

### TODO: getBiomartAttribGroups() is redundant with the above global
### constants. They need to be merged!
getBiomartAttribGroups <- function(id_prefix="ensembl_") {
  attribs <- list()
  attribs[['A1']] <- c(paste0(id_prefix, "transcript_id"),
                       "chromosome_name",
                       "strand",
                       "transcript_start",
                       "transcript_end")

  attribs[['A2']] <- c(paste0(id_prefix, "transcript_id"),
                       "strand",
                       "rank",
                       "exon_chrom_start",
                       "exon_chrom_end")

  attribs[['B']] <- paste0(id_prefix, "exon_id")

  attribs[['C1']] <- c("genomic_coding_start",
                       "genomic_coding_end")

  attribs[['C2']] <- c("5_utr_start",
                       "5_utr_end",
                       "3_utr_start",
                       "3_utr_end")

  attribs[['D']] <- c("cds_start",
                      "cds_end",
                      "cds_length")

  attribs[['G']] <- paste0(id_prefix, "gene_id")
  return(attribs)
}

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
    ## Groups B, C1, C2, and D are optional attributes.
    ## Either group C1 or C2 is required for inferring the CDS (they cannot
    ## be inferred from D). If C1 and C2 are missing, the TxDb object can
    ## still be made but won't have any CDS (no row in the cds table).
    if (.B_ATTRIB %in% attribs$name)
        attrgroups <- paste0(attrgroups, "B")
    if (all(.C1_ATTRIBS %in% attribs$name))
        attrgroups <- paste0(attrgroups, "C1")
    if (all(.C2_ATTRIBS %in% attribs$name))
        attrgroups <- paste0(attrgroups, "C2")
    if (all(.D_ATTRIBS %in% attribs$name))
        attrgroups <- paste0(attrgroups, "D")
    ## Group G: Required attribute.
    if (.G_ATTRIB %in% attribs$name)
        attrgroups <- paste0(attrgroups, "G")
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
        tmp <- paste0(names(tbl2), ":", tbl2, collapse=", ")
        cat("table of attribute groups: ", tmp, "\n", sep="")
    }
    cat("\n")
}

findCompatibleMarts <- function(marts=NULL)
{
    if (is.null(marts))
        marts <- listMarts()
    biomarts <- as.character(marts$biomart)
    versions <- as.character(marts$version)
    for (i in seq_len(nrow(marts)))
        scanMart(biomarts[i], versions[i])
}


### =========================================================================
### findCompatibleMarts() output as of 9/23/2014
### Notes:
###  - Only biomarts containing datasets with at least attribute groups A and
###    G are listed.
###  - Since Ensembl Genomes release 17 (Feb 2013), BioMart access to
###    Ensembl Bacteria is no longer possible. See announcement:
###      http://www.ensembl.info/blog/2013/02/07/ensembl-genomes-release-17/
### -------------------------------------------------------------------------
###
### biomart: ensembl
### version: ENSEMBL GENES 75 (SANGER UK)
### nb of datasets: 66 (hsapiens_gene_ensembl, oanatinus_gene_ensembl,
###                     cporcellus_gene_ensembl, gaculeatus_gene_ensembl, ...)
### table of attribute groups: ABC1C2DG:66
###
### biomart: fungi_mart_22
### version: ENSEMBL FUNGI 22 (EBI UK)
### nb of datasets: 45 (aterreus_eg_gene, mlaricipopulina_eg_gene,
###                     treesei_eg_gene, ...)
### table of attribute groups: ABC1C2G:45
###
### biomart: metazoa_mart_22
### version: ENSEMBL METAZOA 22 (EBI UK)
### nb of datasets: 52 (tcastaneum_eg_gene, dgrimshawi_eg_gene,
###                     rprolixus_eg_gene, ...)
### table of attribute groups: ABC1C2G:52
###
### biomart: plants_mart_22
### version: ENSEMBL PLANTS 22 (EBI UK)
### nb of datasets: 33 (atauschii_eg_gene, obrachyantha_eg_gene,
###                     ppersica_eg_gene, ...)
### table of attribute groups: ABC1C2G:33
###
### biomart: protists_mart_22
### version: ENSEMBL PROTISTS 22 (EBI UK)
### nb of datasets: 29 (pramorum_eg_gene, pvivax_eg_gene,
###                     piwayamai_eg_gene, ...)
### table of attribute groups: ABC1C2G:29
###
### biomart: Breast_mart_69
### version: BCCTB Bioinformatics Portal (UK and Ireland)
### nb of datasets: 2 (breastCancer_expressionStudy, hsapiens_gene_breastCancer)
### table of attribute groups: ABC1C2DG:1, none:1
###
### biomart: Pancreas63
### version: PANCREATIC EXPRESSION DATABASE (BARTS CANCER INSTITUTE UK)
### nb of datasets: 2 (hsapiens_gene_pancreas, hsapiens_cancer_pancreas)
### table of attribute groups: ABC1C2DG:1, none:1
###
### biomart: ENSEMBL_MART_PLANT
### version: GRAMENE 40 ENSEMBL GENES (CSHL/CORNELL US)
### nb of datasets: 33 (atauschii_eg_gene, obrachyantha_eg_gene,
###                     ptrichocarpa_eg_gene, ...)
### table of attribute groups: ABC1C2G:33


### =========================================================================
### findCompatibleMarts() output as of 6/28/2010
###
### (Only biomarts containing datasets with at least attribute groups A and G
### are listed.)
### -------------------------------------------------------------------------
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

