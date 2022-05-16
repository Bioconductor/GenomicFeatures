### =========================================================================
### makeTxDbFromGFF()
### -------------------------------------------------------------------------


.make_Seqinfo_from_chrominfo <- function(chrominfo, circ_seqs=NULL)
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
        is_circ <- GenomeInfoDb:::make_circ_flags_from_circ_seqs(
                                            seqlevels(ans),
                                            circ_seqs)
    } else if (!is.null(circ_seqs)) {
        stop(wmsg("'circ_seqs' should not be specified when 'chrominfo' ",
                  "has an \"is_circular\" column"))
    } else {
        is_circ <- chrominfo$is_circular
    }
    isCircular(ans) <- is_circ
    ans
}

.tidy_seqinfo <- function(gr, circ_seqs=NULL, chrominfo=NULL)
{
    if (is.null(chrominfo)) {
        seqlevels <- seqlevels(gr)
        seqlevels[rankSeqlevels(seqlevels)] <- seqlevels
        seqlevels(gr) <- seqlevels
        isCircular(gr) <- GenomeInfoDb:::make_circ_flags_from_circ_seqs(
                                                   seqlevels(gr),
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
                                taxonomyId=NA, miRBaseBuild=NA, metadata)
{
    message("Prepare the 'metadata' data frame ... ", appendLF=FALSE)
    if (!isSingleStringOrNA(dataSource))
        stop("'dataSource' must be a single string or NA")
    if (!isSingleStringOrNA(organism))
        stop("'organism' must be a single string or NA")
    if (!isSingleStringOrNA(miRBaseBuild))
        stop("'miRBaseBuild' must be a single string or NA")
    if (identical(dataSource, NA)) {
        if (is.character(file)) {
            dataSource <- file
        } else {
            dataSource <- as.character(file)
            if (inherits(file, "connection"))
                dataSource <- showConnections(all=TRUE)[dataSource,
                                                        "description"]
        }
    }
    if (!is.na(taxonomyId)) {
        GenomeInfoDb:::check_tax_id(taxonomyId)
    } else if (!is.na(organism)) {
        taxonomyId <- GenomeInfoDb:::lookup_tax_id_by_organism(organism)
    }
    df <- data.frame(
        name=c("Data source", "Organism", "Taxonomy ID", "miRBase build ID"),
        value=c(dataSource, organism, taxonomyId, miRBaseBuild))
    metadata <- rbind(df, metadata)
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
    if (is(file, "BiocFile")) {
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


#' Make a TxDb object from annotations available as a GFF3 or GTF file
#' 
#' The \code{makeTxDbFromGFF} function allows the user to make a \link{TxDb}
#' object from transcript annotations available as a GFF3 or GTF file.
#' 
#' \code{makeTxDbFromGFF} is a convenience function that feeds data from the
#' parsed file to the \code{\link{makeTxDbFromGRanges}} function.
#' 
#' @param file Input GFF3 or GTF file. Can be a path to a file, or an URL, or a
#' connection object, or a \link[rtracklayer]{GFF3File} or
#' \link[rtracklayer]{GTFFile} object.
#' @param format Format of the input file. Accepted values are: \code{"auto"}
#' (the default) for auto-detection of the format, \code{"gff3"}, or
#' \code{"gtf"}.  Use \code{"gff3"} or \code{"gtf"} only if auto-detection
#' failed.
#' @param dataSource A single string describing the origin of the data file.
#' Please be as specific as possible.
#' @param organism What is the Genus and species of this organism. Please use
#' proper scientific nomenclature for example: "Homo sapiens" or "Canis
#' familiaris" and not "human" or "my fuzzy buddy". If properly written, this
#' information may be used by the software to help you out later.
#' @param taxonomyId By default this value is NA and the organism provided will
#' be used to look up the correct value for this. But you can use this argument
#' to override that and supply your own taxonomy id here (which will be
#' separately validated). Since providing a valid taxonomy id will not require
#' us to look up one based on your organism: this is one way that you can
#' loosen the restrictions about what is and isn't a valid value for the
#' organism.
#' @param circ_seqs A character vector to list out which chromosomes should be
#' marked as circular.
#' @param chrominfo Data frame containing information about the chromosomes.
#' Will be passed to the internal call to \code{\link{makeTxDb}}.  See
#' \code{?\link{makeTxDb}} for more information.  Alternatively, can be a
#' \link[GenomeInfoDb]{Seqinfo} object.
#' @param miRBaseBuild Specify the string for the appropriate build Information
#' from mirbase.db to use for microRNAs. This can be learned by calling
#' \code{supportedMiRBaseBuildValues}. By default, this value will be set to
#' \code{NA}, which will inactivate the \code{microRNAs} accessor.
#' @param metadata A 2-column data frame containing meta information to be
#' included in the \link{TxDb} object. See \code{?\link{makeTxDb}} for more
#' information about the format of \code{metadata}.
#' @param dbxrefTag If not missing, the values in the \code{Dbxref} attribute
#' with the specified tag (like \dQuote{GeneID}) are used for the feature
#' names.
#' @return A \link{TxDb} object.
#' @author M. Carlson and H. Pagès
#' @seealso \itemize{ \item \code{\link{makeTxDbFromGRanges}}, which
#' \code{makeTxDbFromGFF} is based on, for making a \link{TxDb} object from a
#' \link[GenomicRanges]{GRanges} object.
#' 
#' \item The \code{\link[BiocIO]{import}} function in the \pkg{rtracklayer}
#' package (also used by \code{makeTxDbFromGFF} internally).
#' 
#' \item \code{\link{makeTxDbFromUCSC}}, \code{\link{makeTxDbFromBiomart}}, and
#' \code{\link{makeTxDbFromEnsembl}}, for making a \link{TxDb} object from
#' online resources.
#' 
#' \item The \code{\link{supportedMiRBaseBuildValues}} function for listing all
#' the possible values for the \code{miRBaseBuild} argument.
#' 
#' \item The \link{TxDb} class.
#' 
#' \item \code{\link{makeTxDb}} for the low-level function used by the
#' \code{makeTxDbFrom*} functions to make the \link{TxDb} object returned to
#' the user.  }
#' @examples
#' 
#' ## TESTING GFF3
#' gffFile <- system.file("extdata","GFF3_files","a.gff3",package="GenomicFeatures")
#' txdb <- makeTxDbFromGFF(file=gffFile,
#'             dataSource="partial gtf file for Tomatoes for testing",
#'             organism="Solanum lycopersicum")
#' 
#' ## TESTING GTF, this time specifying the chrominfo
#' gtfFile <- system.file("extdata","GTF_files","Aedes_aegypti.partial.gtf",
#'                        package="GenomicFeatures")
#' chrominfo <- data.frame(chrom = c('supercont1.1','supercont1.2'),
#'                         length=c(5220442, 5300000),
#'                         is_circular=c(FALSE, FALSE))
#' metadata <- data.frame(name="Resource URL",
#'                        value=paste0("ftp://ftp.ensemblgenomes.org/pub/metazoa/",
#'                                     "release-13/gtf/aedes_aegypti/"))
#' txdb2 <- makeTxDbFromGFF(file=gtfFile,
#'              chrominfo=chrominfo,
#'              dataSource="ensemblgenomes",
#'              organism="Aedes aegypti",
#'              metadata=metadata)
#' 
#' @export makeTxDbFromGFF
makeTxDbFromGFF <- function(file,
                            format=c("auto", "gff3", "gtf"),
                            dataSource=NA,
                            organism=NA,
                            taxonomyId=NA,
                            circ_seqs=NULL,
                            chrominfo=NULL,
                            miRBaseBuild=NA,
                            metadata=NULL,
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
                                    miRBaseBuild, metadata)

    message("Make the TxDb object ... ", appendLF=FALSE)
    txdb <- makeTxDbFromGRanges(gr, metadata=metadata)
    message("OK")
    txdb
}

