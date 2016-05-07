### =========================================================================
### makeTxDbFromGFF()
### -------------------------------------------------------------------------


.make_Seqinfo_from_chrominfo <- function(chrominfo,
                                         circ_seqs=DEFAULT_CIRC_SEQS)
{
    if (is(chrominfo, "Seqinfo"))
        return(chrominfo)
    if (!is.data.frame(chrominfo))
        stop(wmsg("'chrominfo' must be a data.frame, a Seqinfo object, ",
                  "or NULL"))
    .REQUIRED_COLS <- c("chrom", "length")
    .OPTIONAL_COLS <- "is_circular"
    check_colnames(chrominfo, .REQUIRED_COLS, .OPTIONAL_COLS, "chrominfo")
    ans <- Seqinfo(as.character(chrominfo$chrom), chrominfo$length)
    if (!has_col(chrominfo, "is_circular")) {
        is_circ <- make_circ_flags_from_circ_seqs(seqlevels(ans), circ_seqs)
    } else if (!identical(circ_seqs, DEFAULT_CIRC_SEQS)) {
        stop(wmsg("'circ_seqs' should not be specified when 'chrominfo' ",
                  "has an \"is_circular\" column"))
    } else {
        is_circ <- chrominfo$is_circular
    }
    isCircular(ans) <- is_circ
    ans
}

.tidy_seqinfo <- function(gr, circ_seqs=DEFAULT_CIRC_SEQS, chrominfo=NULL)
{
    if (is.null(chrominfo)) {
        seqlevels <- seqlevels(gr)
        seqlevels[rankSeqlevels(seqlevels)] <- seqlevels
        seqlevels(gr) <- seqlevels
        isCircular(gr) <- make_circ_flags_from_circ_seqs(seqlevels(gr),
                                                         circ_seqs)
        return(gr)
    }
    si <- .make_Seqinfo_from_chrominfo(chrominfo, circ_seqs)
    suppressWarnings(seqinfo(gr) <- merge(seqinfo(gr), si))
    dangling_seqlevels <- setdiff(seqlevelsInUse(gr), seqlevels(si))
    if (length(dangling_seqlevels) != 0L) {
        in1string <- paste0(dangling_seqlevels, collapse=", ")
        stop(wmsg("'chrominfo' must describe at least all the chromosomes ",
                  "of the genomic features imported from the file. ",
                  "Chromosomes missing from 'chrominfo': ", in1string))
    }
    seqlevels(gr) <- seqlevels(si)
    gr
}

.rename_by_dbxrefTag <- function(gr, dbxrefTag) {
    dbxref <- unlist(gr$Dbxref, use.names=FALSE)
    m <- grepl(paste0("^", dbxrefTag, ":"), dbxref)
    gr$Name[S4Vectors:::quick_togroup(gr$Dbxref)[m]] <- dbxref[m]
    gr
}

.prepareGFFMetadata <- function(file, dataSource=NA, organism=NA,
                                taxonomyId=NA, miRBaseBuild=NA)
{
    message("Prepare the 'metadata' data frame ... ", appendLF=FALSE)
    if (!isSingleStringOrNA(dataSource))
        stop("'dataSource' must be a a single string or NA")
    if (!isSingleStringOrNA(organism))
        stop("'organism' must be a a single string or NA")
    if (!isSingleStringOrNA(miRBaseBuild))
        stop("'miRBaseBuild' must be a a single string or NA")
    if (identical(dataSource, NA)) {
        if (is.character(file)) {
            dataSource <- file
        } else {
            dataSource <- showConnections(all=TRUE)[as.character(file),
                                                    "description"]
        }
    }
    if(is.na(taxonomyId)){
        taxonomyId <- GenomeInfoDb:::.taxonomyId(organism)
    }else{
        GenomeInfoDb:::.checkForAValidTaxonomyId(taxonomyId)
    }
    metadata <- data.frame(
                   name=c("Data source",
                          "Organism",
                          "Taxonomy ID",
                          "miRBase build ID"),
                   value=c(dataSource, organism,
                     taxonomyId,
                     miRBaseBuild)
                   )
    message("OK")
    metadata
}

.detect_file_format <- function(file)
{
    if (is(file, "connection"))
        file <- summary(file)$description
    if (isSingleString(file)) {
        file2 <- try(FileForFormat(file), silent=TRUE)
        if (inherits(file2, "try-error"))
            return(tools::file_ext(file))
        file <- file2
    }
    if (is(file, "RTLFile")) {
        if (is(file, "GFF3File"))
            return("gff3")
        if (is(file, "GTFFile"))
            return("gtf")
        desc <- rtracklayer:::resourceDescription(file)
        if (is(file, "CompressedFile")) {
            ## Mimic what import,CompressedFile,missing,missing method does.
            desc <- tools::file_path_sans_ext(desc)
        }
        format <- tools::file_ext(desc)
        return(format)
    }
    stop(wmsg("Invalid 'file'. Must be a path to a file, or an URL, ",
              "or a connection object, or a GFF3File or GTFFile object."))
}

### Based on makeTxDbFromGRanges().
### Timing (T) and memory usage (VIRT) of
###   txdb <- makeTxDbFromGFF("dmel-all-r6.03.gff.gz")
### on rhino04 since BioC 3.0:
###   BioC      T    VIRT  comment
###   ----  -----  ------  -------
###    3.0                 old makeTranscriptDbFromGFF()
###    3.1          > 14g  use new makeTxDbFromGRanges() internally
###    3.2    41s  755212  improvements to the import() step (in rtracklayer)
makeTxDbFromGFF <- function(file,
                            format=c("auto", "gff3", "gtf"),
                            dataSource=NA,
                            organism=NA,
                            taxonomyId=NA,
                            circ_seqs=DEFAULT_CIRC_SEQS,
                            chrominfo=NULL,
                            miRBaseBuild=NA,
                            dbxrefTag)
{
    format <- match.arg(format)
    if (format == "auto") {
        format <- .detect_file_format(file)
        if (!(format %in% c("gff3", "gff", "gtf")))
            stop(wmsg("Cannot detect whether 'file' is a GFF3 or GTF file. ",
                      "Please use the 'format' argument to specify the ",
                      "format (\"gff3\" or \"gtf\")."))
    }
    if (format == "gff3") {
        colnames <- GFF3_COLNAMES
    } else if (format == "gtf") {
        colnames <- GTF_COLNAMES
    } else { # format == "gff"
        ## We don't know a priori if the file is GFF3 or GTF.
        ## TODO: Maybe use sniffGFFVersion() to detect whether the file is
        ## GFF3 or GTF (maybe do this in .detect_file_format()).
        colnames <- union(GFF3_COLNAMES, GTF_COLNAMES)
    }

    message("Import genomic features from the file as a GRanges object ... ",
            appendLF=FALSE)
    gr <- import(file, format=format, colnames=colnames,
                       feature.type=GFF_FEATURE_TYPES)
    gr <- .tidy_seqinfo(gr, circ_seqs, chrominfo)
    if (!missing(dbxrefTag)) {
        gr <- .rename_by_dbxrefTag(gr, dbxrefTag)
    }

    message("OK")

    metadata <- .prepareGFFMetadata(file, dataSource, organism, taxonomyId,
                                    miRBaseBuild)

    message("Make the TxDb object ... ", appendLF=FALSE)
    txdb <- makeTxDbFromGRanges(gr, metadata=metadata)
    message("OK")
    txdb
}

