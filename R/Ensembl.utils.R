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


.ENSEMBL.PUB_FTP_URL <- "ftp://ftp.ensembl.org/pub/"

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

### Return URL of Ensemble Core DB (FTP access).
.Ensembl.getMySQLCoreUrl <- function(dataset, release=NA, url=NA)
{
    if (is.na(url))
        url <- .Ensembl.getFtpUrlToMySQL(release)
    core_dir <- .Ensembl.getMySQLCoreDir(dataset, release=release, url=url)
    paste(url, core_dir, "/", sep="")
}

.Ensembl.getTable <- function(base_url, tablename, col.names)
{
    url <- paste(base_url, tablename, ".txt.gz", sep="")
    destfile <- tempfile()
    download.file(url, destfile, quiet=TRUE)
    data <- read.table(destfile, sep="\t", quote="",
                       col.names=col.names, comment.char="",
                       stringsAsFactors=FALSE)
    unlink(destfile)
    data
}

.Ensembl.getTable.seq_region <- function(core_url, with.coord_system=FALSE)
{
    ## Get 'seq_region' table.
    COLNAMES <- c("seq_region_id", "name", "coord_system_id", "length")
    ans <- .Ensembl.getTable(core_url, "seq_region", COLNAMES)
    if (!with.coord_system)
        return(ans)
    ## Get 'coord_system' table.
    COLNAMES <- c("coord_system_id", "species_id", "name",
                  "version", "rank", "attrib")
    COLNAMES[3:6] <- paste("coord_system_", COLNAMES[3:6], sep="")
    coord_system <- .Ensembl.getTable(core_url, "coord_system", COLNAMES)
    left2right <- match(ans[["coord_system_id"]],
                        coord_system[["coord_system_id"]])
    ans_right <- coord_system[left2right, , drop=FALSE]
    rownames(ans_right) <- NULL
    ans <- ans[-match("coord_system_id", colnames(ans))]
    cbind(ans, ans_right)
}

.Ensembl.fetchAttribTypeIdForTopLevelSequence <- function(core_url)
{
    ## Get 'attrib_type' table.
    COLNAMES <- c("attrib_type_id", "code", "name", "description")
    attrib_type <- .Ensembl.getTable(core_url, "attrib_type", COLNAMES)
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
    COLNAMES <- c("seq_region_id", "attrib_type_id", "value")
    seq_region_attrib <- .Ensembl.getTable(core_url,
                             "seq_region_attrib", COLNAMES)
    seq_region_attrib$seq_region_id[seq_region_attrib$attrib_type_id == id0]
}

### Fetch sequence names and lengths from the 'seq_region' table.
### Typical use:
###   core_url <- .Ensembl.getMySQLCoreUrl("hsapiens_gene_ensembl")
###   extra_seqnames <- c("GL000217.1", "NC_012920", "HG79_PATCH")
###   .Ensembl.fetchChromLengthsFromCoreUrl(core_url,
###                                         extra_seqnames=extra_seqnames)
.Ensembl.fetchChromLengthsFromCoreUrl <- function(core_url, extra_seqnames=NULL)
{
    seq_region <- .Ensembl.getTable.seq_region(core_url, with.coord_system=TRUE)

    ## 1st filtering: Keep only "default_version" sequences.
    i1 <- grep("default_version", seq_region$coord_system_attrib, fixed=TRUE)
    j1 <- c("seq_region_id", "name", "length",
            "coord_system_name", "coord_system_rank")
    ans <- seq_region[i1, j1, drop=FALSE]

    ## 2nd filtering: Keep only "toplevel" sequences that are not LRGs +
    ## extra sequences.
    ids <- .Ensembl.fetchTopLevelSequenceIds(core_url)
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

    ## Ordering: First by rank, then by name (main chromosomes first, extra
    ## sequences last).
    names_as_int <- suppressWarnings(as.integer(ans$name))
    nb_ints <- sum(!is.na(names_as_int))
    names_as_int[match("X", ans$name)] <- nb_ints + 1L
    names_as_int[match("Y", ans$name)] <- nb_ints + 2L
    names_as_int[match("M", ans$name)] <- nb_ints + 3L
    names_as_int[match("MT", ans$name)] <- nb_ints + 4L
    names_as_int[match(extra_seqnames, ans$name)] <- .Machine$integer.max
    prev_locale <- Sys.getlocale("LC_COLLATE")
    Sys.setlocale("LC_COLLATE", "C")
    oo <- order(ans$coord_system_rank, names_as_int, ans$name)
    Sys.setlocale("LC_COLLATE", prev_locale)
    ans <- ans[oo, , drop=FALSE]

    ## 3rd filtering: There can be more than one row per sequence name, but
    ## the pair (name, coord_system_rank) should be unique. We disambiguate
    ## by keeping rows with the lowest coord_system_rank. This is
    ## straightforward because the rows are already ordered from lowest to
    ## highest ranks.
    i3 <- !duplicated(ans$name)
    j3 <- c("name", "length")
    ans <- ans[i3, j3, drop=FALSE]

    ## Final tidying.
    row.names(ans) <- NULL
    ans
}

fetchChromLengthsFromEnsembl <- function(dataset, release=NA,
                                         extra_seqnames=NULL)
{
    core_url <- .Ensembl.getMySQLCoreUrl(dataset, release=release)
    .Ensembl.fetchChromLengthsFromCoreUrl(core_url,
                                          extra_seqnames=extra_seqnames)
}

