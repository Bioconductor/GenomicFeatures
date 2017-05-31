### =========================================================================
### Utilities for direct query of the Ensembl FTP server
### -------------------------------------------------------------------------
###
### The utililities in this file use RCurl or just utils::download.file() for
### getting stuff directly from the Ensembl FTP server. They can access stuff
### that is not available thru biomaRt like for example the lengths of the
### sequences in the reference genome associated with a particular dataset
### and Ensembl release (e.g. for dataset "hsapiens_gene_ensembl" and Ensembl
### release "64").
### Note that querying the Ensembl MySQL server (via RMySQL) would probably
### be a better way to access this stuff but that would mean one more
### dependency for the GenomicFeatures package. With some potential
### complications like: (a) no RMySQL Windows binary on CRAN, and (b) depending
### on RMySQL *and* RSQLite has its own pitfalls.
###
### Ensembl Core Schema Documentation:
###   http://www.ensembl.org/info/docs/api/core/core_schema.html
### The full schema:
###   ftp://ftp.ensembl.org/pub/ensembl/sql/table.sql
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Very low-level internal utils used here and in other files
###

### Uses RCurl to access and list the content of an FTP dir.
ls_ftp_url <- function(url, subdirs.only=FALSE)
{
    doc <- getURL(url)
    listing <- strsplit(doc, "\n", fixed=TRUE)[[1L]]
    if (subdirs.only)
        listing <- listing[substr(listing, 1L, 1L) == "d"]
    ## Keep field no. 8 only
    pattern <- paste(c("^", rep.int("[^[:space:]]+[[:space:]]+", 8L)),
                     collapse="")
    listing <- sub(pattern, "", listing)
    sub("[[:space:]].*$", "", listing)
}

.ENSEMBL.PUB_FTP_URL <- "ftp://ftp.ensembl.org/pub/"
.ENSEMBLGRCh37.PUB_FTP_URL <- "ftp://ftp.ensembl.org/pub/grch37/"

ftp_url_to_Ensembl_mysql <- function(release=NA, use.grch37=FALSE)
{
    if (is.na(release)) {
        if (use.grch37)
            pub_subdir <- "current/mysql"
        else 
            pub_subdir <- "current_mysql"
    } else {
        pub_subdir <- paste0("release-", release, "/mysql")
    }
    if (use.grch37)
        pub_ftp_url <- .ENSEMBLGRCh37.PUB_FTP_URL
    else
        pub_ftp_url <- .ENSEMBL.PUB_FTP_URL
    paste0(pub_ftp_url, pub_subdir, "/")
}

ftp_url_to_Ensembl_gtf <- function(release=NA)
{
    if (is.na(release))
        pub_subdir <- "current_gtf"
    else
        pub_subdir <- paste0("release-", release, "/gtf")
    paste0(.ENSEMBL.PUB_FTP_URL, pub_subdir, "/")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .Ensembl_getMySQLCoreUrl()
###

.Ensembl_listMySQLCoreDirs <- function(release=NA, url=NA, use.grch37=FALSE)
{
    if (is.na(url))
        url <- ftp_url_to_Ensembl_mysql(release, use.grch37=use.grch37)
    core_dirs <- ls_ftp_url(url, subdirs.only=TRUE)
    pattern <- "_core_"
    if (!is.na(release))
        pattern <- paste0(pattern, release, "_")
    core_dirs[grep(pattern, core_dirs, fixed=TRUE)]
}

.Ensembl_getMySQLCoreDir <- function(dataset, release=NA, url=NA,
                                     use.grch37=FALSE)
{
    if (is.na(url))
        url <- ftp_url_to_Ensembl_mysql(release, use.grch37=use.grch37)
    core_dirs <- .Ensembl_listMySQLCoreDirs(release=release, url=url,
                                            use.grch37=use.grch37)
    trimmed_core_dirs <- sub("_core_.*$", "", core_dirs)
    shortnames <- sub("^(.)[^_]*_", "\\1", trimmed_core_dirs)
    if (dataset == "mfuro_gene_ensembl") {
        shortname0 <- "mputorius_furo"
    } else {
        shortname0 <- strsplit(dataset, "_", fixed=TRUE)[[1L]][1L]
    }
    core_dir <- core_dirs[shortnames == shortname0]
    if (length(core_dir) != 1L)
        stop("found 0 or more than 1 subdir for \"", dataset,
             "\" dataset at ", url)
    core_dir
}

### Return URL of Ensemble Core DB (FTP access).
.Ensembl_getMySQLCoreUrl <- function(dataset, release=NA, url=NA,
                                     use.grch37=FALSE)
{
    if (is.na(url))
        url <- ftp_url_to_Ensembl_mysql(release, use.grch37=use.grch37)
    core_dir <- .Ensembl_getMySQLCoreDir(dataset, release=release, url=url,
                                         use.grch37=use.grch37)
    paste0(url, core_dir, "/")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .Ensembl_fetchChromLengthsFromCoreUrl()
###

.Ensembl_getTable <- function(base_url, tablename, col.names, nrows=-1)
{
    url <- paste0(base_url, tablename, ".txt.gz")
    destfile <- tempfile()
    download.file(url, destfile, quiet=TRUE)
    data <- read.table(destfile, sep="\t", quote="",
                       col.names=col.names, nrows=nrows, comment.char="",
                       stringsAsFactors=FALSE)
    unlink(destfile)
    data
}

.Ensembl_getTable_seq_region <- function(core_url, with.coord_system=FALSE)
{
    ## Get 'seq_region' table.
    COLNAMES <- c("seq_region_id", "name", "coord_system_id", "length")
    ans <- .Ensembl_getTable(core_url, "seq_region", COLNAMES)
    if (!with.coord_system)
        return(ans)
    ## Get 'coord_system' table.
    COLNAMES <- c("coord_system_id", "species_id", "name",
                  "version", "rank", "attrib")
    COLNAMES[3:6] <- paste0("coord_system_", COLNAMES[3:6])
    coord_system <- .Ensembl_getTable(core_url, "coord_system", COLNAMES)
    left2right <- match(ans[["coord_system_id"]],
                        coord_system[["coord_system_id"]])
    ans_right <- S4Vectors:::extract_data_frame_rows(coord_system, left2right)
    ans <- ans[-match("coord_system_id", colnames(ans))]
    cbind(ans, ans_right)
}

.Ensembl_fetchAttribTypeIdForTopLevelSequence <- function(core_url)
{
    ## Get the first 50 rows of 'attrib_type' table.
    ## The reason we don't read in the entire table is because some rows at
    ## the bottom of the table (rows with attrib_type_id >= 416, these rows
    ## were added in Ensembl 75) contain embedded EOL characters that break
    ## read.table(). Since we only need to retrieve the attrib_type_id
    ## associated with the "toplevel" code, and since this code is generally
    ## found at the beginning of the file (AFAIK "toplevel" has always been
    ## assigned the attrib_type_id of 6 and the lines in the file seem to
    ## always be ordered by attrib_type_id), reading in the first 50 rows
    ## should be way enough to get what we need.
    COLNAMES <- c("attrib_type_id", "code", "name", "description")
    attrib_type <- .Ensembl_getTable(core_url, "attrib_type", COLNAMES,
                                     nrows=50)
    i <- which(attrib_type$code == "toplevel")
    if (length(i) != 1L)
        stop("Ensembl data anomaly: \"toplevel\" attrib found 0 or more ",
             "than once in attrib_type table at ", core_url)
    attrib_type$attrib_type_id[i]
}

.Ensembl_fetchTopLevelSequenceIds <- function(core_url)
{
    id0 <- .Ensembl_fetchAttribTypeIdForTopLevelSequence(core_url)
    ## Get 'seq_region_attrib' table.
    COLNAMES <- c("seq_region_id", "attrib_type_id", "value")
    seq_region_attrib <- .Ensembl_getTable(core_url,
                             "seq_region_attrib", COLNAMES)
    seq_region_attrib$seq_region_id[seq_region_attrib$attrib_type_id == id0]
}

### Fetch sequence names and lengths from the 'seq_region' table.
### Typical use:
###   core_url <- .Ensembl_getMySQLCoreUrl("hsapiens_gene_ensembl")
###   extra_seqnames <- c("GL000217.1", "NC_012920", "HG79_PATCH")
###   .Ensembl_fetchChromLengthsFromCoreUrl(core_url,
###                                         extra_seqnames=extra_seqnames)
.Ensembl_fetchChromLengthsFromCoreUrl <- function(core_url, extra_seqnames=NULL)
{
    seq_region <- .Ensembl_getTable_seq_region(core_url, with.coord_system=TRUE)

    ## 1st filtering: Keep only "default_version" sequences.
    i1 <- grep("default_version", seq_region$coord_system_attrib, fixed=TRUE)
    j1 <- c("seq_region_id", "name", "length",
            "coord_system_name", "coord_system_rank")
    ans <- seq_region[i1, j1, drop=FALSE]

    ## 2nd filtering: Keep only "toplevel" sequences that are not LRGs +
    ## extra sequences.
    ids <- .Ensembl_fetchTopLevelSequenceIds(core_url)
    i2 <- ans$seq_region_id %in% ids & ans$coord_system_name != "lrg"
    if (!is.null(extra_seqnames)) {
        extra_seqnames <- unique(extra_seqnames)
        if (!all(extra_seqnames %in% ans$name))
            stop("failed to fetch all chromosome lengths")
        extra_seqnames <- setdiff(extra_seqnames, ans$name[i2])
        ## Add extra sequences to the index.
        i2 <- i2 | (ans$name %in% extra_seqnames)
    }
    j2 <- c("name", "length", "coord_system_rank")
    ans <- ans[i2, j2, drop=FALSE]

    ## Ordering: First by rank, then by name.
    names_as_int <- rankSeqlevels(ans$name)
    oo <- order(ans$coord_system_rank, names_as_int)
    ans <- S4Vectors:::extract_data_frame_rows(ans, oo)

    ## 3rd filtering: There can be more than one row per sequence name, but
    ## the pair (name, coord_system_rank) should be unique. We disambiguate
    ## by keeping rows with the lowest coord_system_rank. This is
    ## straightforward because the rows are already ordered from lowest to
    ## highest ranks.
    i3 <- !duplicated(ans$name)
    j3 <- c("name", "length")
    ans <- ans[i3, j3, drop=FALSE]

    ## Final tidying.
    rownames(ans) <- NULL
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_organism_from_Ensembl_Mart_dataset()
###

get_organism_from_Ensembl_Mart_dataset <- function(dataset)
{
    core_dir <- .Ensembl_getMySQLCoreDir(dataset)
    organism <- sub("_core.*", "", core_dir)
    organism <- sub("_", " ", organism)
    substr(organism, 1L, 1L) <- toupper(substr(organism, 1L, 1L))
    organism
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetchChromLengthsFromEnsembl() and fetchChromLengthsFromEnsemblPlants()
###

fetchChromLengthsFromEnsembl <- function(dataset, release=NA, use.grch37=FALSE,
                                         extra_seqnames=NULL)
{
    core_url <- .Ensembl_getMySQLCoreUrl(dataset, release=release,
                                         use.grch37=use.grch37)
    .Ensembl_fetchChromLengthsFromCoreUrl(core_url,
                                          extra_seqnames=extra_seqnames)
}

.ENSEMBL_PLANTS.CURRENT_MYSQL_URL <- "ftp://ftp.ensemblgenomes.org/pub/plants/current/mysql/"

fetchChromLengthsFromEnsemblPlants <- function(dataset,
                                               extra_seqnames=NULL)
{
    core_url <- .Ensembl_getMySQLCoreUrl(dataset,
                                         url=.ENSEMBL_PLANTS.CURRENT_MYSQL_URL)
    .Ensembl_fetchChromLengthsFromCoreUrl(core_url,
                                          extra_seqnames=extra_seqnames)
}

