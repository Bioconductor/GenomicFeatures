### =========================================================================
### Utilities for direct query of the Ensembl FTP server
### -------------------------------------------------------------------------
###
### The utililities in this file use utils::download.file() to get stuff
### directly from the Ensembl FTP server. They can access stuff that is not
### available thru biomaRt like for example the lengths of the sequences in
### the reference genome associated with a particular dataset and Ensembl
### release (e.g. for dataset "hsapiens_gene_ensembl" and Ensembl release "64").
### Note that querying the Ensembl MySQL server (via RMariaDB) would probably
### be a better way to access this stuff.
###
### Ensembl Core Schema Documentation:
###   http://www.ensembl.org/info/docs/api/core/core_schema.html
### Full schema:
###   ftp://ftp.ensembl.org/pub/ensembl/sql/table.sql
###

.ENSEMBL.PUB_FTP_URL <- "ftp://ftp.ensembl.org/pub/"
.ENSEMBLGRCh37.PUB_FTP_URL <- "ftp://ftp.ensembl.org/pub/grch37/"
.ENSEMBLGENOMES.PUB_FTP_URL <- "ftp://ftp.ensemblgenomes.org/pub/"


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Very low-level internal utils used here and in other files
###

.normarg_release <- function(release)
{
    error_msg <- "'release' must be a valid Ensembl release number or NA"
    if (!(is.atomic(release) && length(release) == 1L))
        stop(wmsg(error_msg))
    if (!(is.integer(release) || is.na(release))) {
        if (!(is.numeric(release) || is.character(release)))
            stop(wmsg(error_msg))
        release2 <- suppressWarnings(as.integer(release))
        if (is.na(release2) || release2 != release)
            stop(wmsg(error_msg))
        release <- release2
    }
    release
}

### 'kingdom' must be NA or one of the EnsemblGenomes marts i.e. "bacteria",
### "fungi", "metazoa", "plants", or "protists".
ftp_url_to_Ensembl_mysql <- function(release=NA, use.grch37=FALSE, kingdom=NA)
{
    release <- .normarg_release(release)
    if (is.na(kingdom)) {
        if (is.na(release)) {
            if (use.grch37) {
                pub_subdir <- "current/mysql"
            } else {
                pub_subdir <- "current_mysql"
            }
        } else {
            pub_subdir <- paste0("release-", release, "/mysql")
        }
        if (use.grch37)
            pub_ftp_url <- .ENSEMBLGRCh37.PUB_FTP_URL
        else
            pub_ftp_url <- .ENSEMBL.PUB_FTP_URL
    } else {
        pub_ftp_url <- paste0(.ENSEMBLGENOMES.PUB_FTP_URL, kingdom, "/")
        if (is.na(release)) {
            pub_subdir <- "current"
        } else {
            pub_subdir <- paste0("release-", release)
        }
        pub_subdir <- paste0(pub_subdir, "/mysql")
    }
    paste0(pub_ftp_url, pub_subdir, "/")
}

ftp_url_to_Ensembl_gtf <- function(release=NA)
{
    release <- .normarg_release(release)
    if (is.na(release))
        pub_subdir <- "current_gtf"
    else
        pub_subdir <- paste0("release-", release, "/gtf")
    paste0(.ENSEMBL.PUB_FTP_URL, pub_subdir, "/")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .Ensembl_getMySQLCoreUrl()
###

Ensembl_listMySQLCoreDirs <- function(mysql_url, release=NA, recache=FALSE)
{
    release <- .normarg_release(release)
    core_dirs <- list_ftp_dir(mysql_url, subdirs.only=TRUE, recache=recache)
    pattern <- "_core_"
    if (!is.na(release))
        pattern <- paste0(pattern, release, "_")
    core_dirs[grep(pattern, core_dirs, fixed=TRUE)]
}

.Ensembl_getMySQLCoreDir <- function(dataset, mysql_url, release=NA)
{
    core_dirs <- Ensembl_listMySQLCoreDirs(mysql_url, release=release)
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
             "\" dataset at ", mysql_url)
    core_dir
}

### Return URL of Ensemble Core DB (FTP access).
.Ensembl_getMySQLCoreUrl <- function(dataset, mysql_url, release=NA)
{
    core_dir <- .Ensembl_getMySQLCoreDir(dataset, mysql_url, release=release)
    paste0(mysql_url, core_dir, "/")
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

extract_chromlengths_from_seq_region <- function(seq_region,
                                                 top_level_ids,
                                                 extra_seqnames=NULL,
                                                 seq_region_ids=NULL)
{
    ## 1st filtering: Keep only "default_version" sequences.
    keep_me <- grepl("default_version", seq_region$coord_system_attrib,
                    fixed=TRUE)
    if (!is.null(seq_region_ids))
        keep_me <- keep_me | (seq_region$seq_region_id %in% seq_region_ids)
    i1 <- which(keep_me)
    j1 <- c("seq_region_id", "name", "length",
            "coord_system_name", "coord_system_rank")
    ans <- seq_region[i1, j1, drop=FALSE]

    ## 2nd filtering: Keep only "toplevel" sequences that are not LRGs +
    ## extra sequences.
    keep_me <- ans$seq_region_id %in% top_level_ids &
               ans$coord_system_name != "lrg"
    if (!is.null(extra_seqnames)) {
        extra_seqnames <- unique(extra_seqnames)
        if (!all(extra_seqnames %in% ans$name))
            stop("failed to fetch all chromosome lengths")
        extra_seqnames <- setdiff(extra_seqnames, ans$name[keep_me])
        ## Add extra sequences to the index.
        keep_me <- keep_me | (ans$name %in% extra_seqnames)
    }
    if (!is.null(seq_region_ids))
        keep_me <- keep_me | (ans$seq_region_id %in% seq_region_ids)
    i2 <- which(keep_me)
    j2 <- c("seq_region_id", "name", "length", "coord_system_rank")
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
    keep_me <- !duplicated(ans$name)
    if (!is.null(seq_region_ids))
        keep_me <- keep_me | (ans$seq_region_id %in% seq_region_ids)
    i3 <- which(keep_me)
    j3 <- c("seq_region_id", "name", "length")
    ans <- ans[i3, j3, drop=FALSE]

    ## Final tidying.
    rownames(ans) <- NULL
    ans
}

### Fetch sequence names and lengths from the 'seq_region' table.
### Typical use:
###   mysql_url <- ftp_url_to_Ensembl_mysql()
###   core_url <- .Ensembl_getMySQLCoreUrl("hsapiens_gene_ensembl", mysql_url)
###   extra_seqnames <- c("GL000217.1", "NC_012920", "HG79_PATCH")
###   .Ensembl_fetchChromLengthsFromCoreUrl(core_url,
###                                         extra_seqnames=extra_seqnames)
.Ensembl_fetchChromLengthsFromCoreUrl <- function(core_url, extra_seqnames=NULL)
{
    seq_region <- .Ensembl_getTable_seq_region(core_url, with.coord_system=TRUE)

    top_level_ids <- .Ensembl_fetchTopLevelSequenceIds(core_url)
    ans <- extract_chromlengths_from_seq_region(seq_region,
                                                top_level_ids,
                                                extra_seqnames=extra_seqnames)
    ans[-1L]  # drop "seq_region_id" col
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### get_organism_from_Ensembl_Mart_dataset()
###

get_organism_from_Ensembl_Mart_dataset <- function(dataset, release=NA,
                                                   use.grch37=FALSE,
                                                   kingdom=NA)
{
    mysql_url <- ftp_url_to_Ensembl_mysql(release, use.grch37, kingdom)
    core_dir <- .Ensembl_getMySQLCoreDir(dataset, mysql_url, release=release)
    organism <- sub("_core.*", "", core_dir)
    organism <- sub("_", " ", organism)
    substr(organism, 1L, 1L) <- toupper(substr(organism, 1L, 1L))
    organism
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### fetchChromLengthsFromEnsembl()
###

### 'kingdom' must be NA or one of the EnsemblGenomes marts i.e. "bacteria",
### "fungi", "metazoa", "plants", or "protists".
fetchChromLengthsFromEnsembl <- function(dataset, release=NA,
                                         use.grch37=FALSE, kingdom=NA,
                                         extra_seqnames=NULL)
{
    mysql_url <- ftp_url_to_Ensembl_mysql(release, use.grch37, kingdom)
    core_url <- .Ensembl_getMySQLCoreUrl(dataset, mysql_url, release=release)
    .Ensembl_fetchChromLengthsFromCoreUrl(core_url,
                                          extra_seqnames=extra_seqnames)
}

