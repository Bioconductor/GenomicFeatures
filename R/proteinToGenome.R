### =========================================================================
### proteinToGenome()
### -------------------------------------------------------------------------


#' Map protein-relative coordinates to genomic coordinates
#'
#' \code{proteinToGenome} is a generic function for mapping ranges of
#' protein-relative positions to the genome.
#'
#' NOTE: This man page is for the \code{proteinToGenome} S4 generic function
#' and methods defined in the \pkg{GenomicFeatures} package, which are
#' (loosely) modeled on the \code{\link[ensembldb]{proteinToGenome}} function
#' from the \pkg{ensembldb} package.  See
#' \code{?ensembldb::\link[ensembldb]{proteinToGenome}} for the latter.
#'
#' The \code{proteinToGenome()} method for \link[GenomicRanges]{GRangesList}
#' objects is the workhorse behind the default method. Note that the latter is
#' a thin wrapper around the former, which simply does the following:
#' \enumerate{ \item Use \code{\link{cdsBy}()} to extract the CDS from
#' \code{txdb}.  The CDS are returned in a \link[GenomicRanges]{GRangesList}
#' object with transcript names on it.  \item Call \code{proteinToGenome()} on
#' \code{x} and the \link[GenomicRanges]{GRangesList} object returned by
#' \code{\link{cdsBy}()}.  }
#'
#' @aliases proteinToGenome proteinToGenome,GRangesList-method
#' proteinToGenome,ANY-method
#' @param x A named \link[IRanges]{IRanges} object (or derivative) containing
#' ranges of \emph{protein-relative positions} (protein-relative positions are
#' positions relative to a protein sequence).
#'
#' The names on \code{x} must be transcript names present in \code{txdb}.  More
#' precisely, for the default \code{proteinToGenome()} method, \code{names(x)}
#' must be a subset of: \preformatted{ mcols(transcripts(txdb,
#' columns="tx_name"))$tx_name } And for the method for
#' \link[GenomicRanges]{GRangesList} objects, \code{names(x)} must be a subset
#' of: \preformatted{ names(txdb) }
#' @param txdb For the default \code{proteinToGenome()} method: A \link{TxDb}
#' object or any object that supports \code{\link{transcripts}()} and
#' \code{\link{cdsBy}()} (e.g. an \link[ensembldb]{EnsDb} object from the
#' \pkg{ensembldb} package).
#'
#' For the method for \link[GenomicRanges]{GRangesList} objects: A named
#' \link[GenomicRanges]{GRangesList} object (or derivative) where each list
#' element is a \link[GenomicRanges]{GRanges} object representing a CDS (the
#' ranges in the \link[GenomicRanges]{GRanges} object must represent the CDS
#' parts ordered by ascending exon rank).
#' @param ...  Further arguments to be passed to specific methods.
#' @return A named \link[GenomicRanges]{GRangesList} object \emph{parallel} to
#' \code{x} (the transcript names on \code{x} are propagated).  The i-th list
#' element in the returned object is the result of mapping the range of
#' protein-relative positions \code{x[i]} to the genome.
#'
#' Note that a given range in \code{x} can only be mapped to the genome if the
#' name on it is the name of a \emph{coding} transcript. If it's not (i.e. if
#' it's the name of a \emph{non-coding} transcript), then an empty
#' \link[GenomicRanges]{GRanges} object is placed in the returned object to
#' indicate the impossible mapping, and a warning is issued.
#'
#' Otherwise, if a given range in \code{x} can be mapped to the genome, then
#' the result of the mapping is represented by a non-empty
#' \link[GenomicRanges]{GRanges} object.  Note that this object represents the
#' original CDS associated to \code{x}, trimmed on its 5' end or 3' end, or on
#' both.  Furthermore, this object will have the same metadata columns as the
#' \link[GenomicRanges]{GRanges} object representing the original CDS, plus the
#' 2 following ones: \itemize{ \item \code{protein_start}: The protein-relative
#' start of the mapping.  \item \code{protein_end}: The protein-relative end of
#' the mapping.  }
#' @note Unlike \code{ensembldb::\link[ensembldb]{proteinToGenome}()} which can
#' work either with Ensembl protein IDs or Ensembl transcript IDs on \code{x},
#' the default \code{proteinToGenome()} method described above only accepts
#' \emph{transcript names} on \code{x}.
#'
#' This means that, if the user is in possession of protein IDs, they must
#' first replace them with the corresponding transcript IDs (referred to as
#' \emph{transcript names} in the context of \link{TxDb} objects).  How to do
#' this exactly depends on the origin of those IDs (UCSC, Ensembl, GTF/GFF3
#' file, FlyBase, etc...)
#' @author H. Pagès, and the authors of \code{ensembldb::proteinToGenome()} for
#' the inspiration and design.
#' @seealso \itemize{ \item The \code{\link[ensembldb]{proteinToGenome}}
#' function in the \pkg{ensembldb} package, which the \code{proteinToGenome()}
#' generic and methods documented in this man page are (loosely) modeled on.
#'
#' \item \link{TxDb} objects.
#'
#' \item \link[ensembldb]{EnsDb} objects (\link{TxDb}-like objects) in the
#' \pkg{ensembldb} package.
#'
#' \item \code{\link{transcripts}} for extracting transcripts from a
#' \link{TxDb}-like object.
#'
#' \item \code{\link{cdsBy}} for extracting CDS from a \link{TxDb}-like object.
#'
#' \item \link[IRanges]{IRanges} objects in the \pkg{IRanges} package.
#'
#' \item \link[GenomicRanges]{GRanges} and \link[GenomicRanges]{GRangesList}
#' objects in the \pkg{GenomicRanges} package.  }
#' @keywords methods utilities
#' @examples
#'
#' ## ---------------------------------------------------------------------
#' ## USING TOY CDS
#' ## ---------------------------------------------------------------------
#' CDS1 <- GRanges(c("chrX:11-60:+", "chrX:101-125:+"))
#' CDS2 <- GRanges(c("chrY:201-230:-", "chrY:101-125:-", "chrY:11-60:-"))
#' cds_by_tx <- GRangesList(TX1=CDS1, TX2=CDS2)
#' cds_by_tx
#'
#' x1 <- IRanges(start=8, end=20, names="TX1")
#' proteinToGenome(x1, cds_by_tx)
#'
#' x2 <- IRanges(start=c(1, 18), end=c(25, 20), names=c("TX1", "TX1"))
#' x2
#' proteinToGenome(x2, cds_by_tx)
#'
#' x3 <- IRanges(start=8, end=15, names="TX2")
#' proteinToGenome(x3, cds_by_tx)
#'
#' x4 <- c(x3, x2)
#' x4
#' proteinToGenome(x4, cds_by_tx)
#'
#' ## ---------------------------------------------------------------------
#' ## USING A TxDb OBJECT
#' ## ---------------------------------------------------------------------
#' library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
#' txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
#'
#' ## The first transcript (FBtr0309810) is non-coding:
#' x <- IRanges(c(FBtr0309810="11-55", FBtr0306539="90-300"))
#' res <- proteinToGenome(x, txdb)
#' res
#'
#' @export
setGeneric("proteinToGenome", signature="txdb",
    function(x, txdb, ...) standardGeneric("proteinToGenome")
)
### Dispatch is on 2nd argument!


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Make fancy error or warning messages
###

.make_bad_names_msg <- function(x_names, bad_idx, what="invalid name",
                                max.show=5L)
{
    nbad <- length(bad_idx)
    if (max.show == 0L) {
        msg <- c("the names on 'x' contain ", nbad, " ", what)
        if (nbad != 1L)
            msg <- c(msg, "s")
        return(paste(msg, collapse=""))
    }
    if (nbad == 1L) {
        msg <- c("The names on 'x' contain ", what)
    } else {
        msg <- c("The names on 'x' contain ", nbad, " ", what, "s")
        if (nbad > max.show) {
            if (max.show == 1L) {
                msg <- c(msg, " (showing the first one only)")
            } else {
                msg <- c(msg, " (showing the first ", max.show, " only)")
            }
            bad_idx <- head(bad_idx, n=max.show)
        }
    }
    bad_names <- x_names[bad_idx]
    paste0(paste(msg, collapse=""), ": ", paste(bad_names, collapse=", "))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### proteinToGenome() method for GRangesList objects
###

.make_protein_ranges_from_cumwidths <- function(cumwidths)
{
    len <- length(cumwidths)
    protein_start <- c(1L, cumwidths[-len] %/% 3L + 1L)
    protein_end <- (2L + cumwidths) %/% 3L
    IRanges(protein_start, protein_end)
}

### Trims first range in 'gr' on its 5' side by 'trim1' nucleotides,
### and last range in 'gr' on its 3' side by 'trim2' nucleotides.
### Other than that, everything else is preserved (length, names, metadata
### columns).
.trim_first_and_last_ranges <- function(gr, trim1=0L, trim2=0L)
{
    stopifnot(is(gr, "GRanges"),
              isSingleInteger(trim1),
              isSingleInteger(trim2))

    gr_ranges <- ranges(gr)
    gr_len <- length(gr_ranges)
    stopifnot(gr_len >= 1L)

    gr_start <- start(gr_ranges)
    gr_end <- end(gr_ranges)
    gr_strand <- S4Vectors:::decodeRle(strand(gr))

    ## Trim first range.
    strand1 <- gr_strand[[1L]]
    if (strand1 == "+") {
        ## Trim on the left.
        gr_start[[1L]] <- gr_start[[1L]] + trim1
    } else {
        ## Trim on the right.
        gr_end[[1L]] <- gr_end[[1L]] - trim1
    }

    ## Trim last range.
    strand2 <- gr_strand[[gr_len]]
    if (strand2 == "+") {
        ## Trim on the right.
        gr_end[[gr_len]] <- gr_end[[gr_len]] - trim2
    } else {
        ## Trim on the left.
        gr_start[[gr_len]] <- gr_start[[gr_len]] + trim2
    }

    if (any(gr_start > gr_end + 1L))
        stop(wmsg("invalid trimming"))

    ranges(gr) <- update_ranges(gr_ranges, start=gr_start, end=gr_end)
    gr
}

### 'cds_parts' must be a GRanges object representing the CDS of a given
### transcript/protein.
### 'protein_start' and 'protein_end' must be protein-relative coordinates
### i.e. coordinates (counted in Amino Acids) relative to the protein
### associated with 'cds_parts'.
### Returns a GRanges object.
.map_protein_to_cds_parts <- function(protein_start, protein_end, cds_parts)
{
    stopifnot(isSingleNumber(protein_start),
              isSingleNumber(protein_end),
              protein_start <= protein_end,
              is(cds_parts, "GRanges"))
    nparts <- length(cds_parts)
    cds_widths <- width(cds_parts)
    stopifnot(nparts >= 1L, all(cds_widths >= 1L))
    protein_start <- as.integer(protein_start)
    protein_end <- as.integer(protein_end)
    cds_cumwidths <- cumsum(cds_widths)

    ## Add metadata columns 'protein_start' and 'protein_end' to 'cds_parts'.
    protein_ranges <- .make_protein_ranges_from_cumwidths(cds_cumwidths)
    protein_ranges <- DataFrame(protein_start=start(protein_ranges),
                                protein_end=end(protein_ranges))
    mcols(cds_parts) <- cbind(mcols(cds_parts), protein_ranges)

    ## Translate protein-relative coordinates into 0-based CDS-relative
    ## coordinates.
    protein_start0 <- 3L * (protein_start - 1L)
    protein_end0 <- 3L * protein_end - 1L

    ## Find CDS parts touched by 'protein_start' and 'protein_end'.
    idx <- 1L + findInterval(c(protein_start0, protein_end0), cds_cumwidths)
    idx1 <- idx[[1L]]
    idx2 <- idx[[2L]]
    if (idx2 > nparts)
        idx2 <- nparts

    ## Extract all CDS parts touched by the [protein_start,protein_end] range.
    ans <- cds_parts[idx1:idx2]

    ## Trim first and last ranges in 'ans' (trimming should **always**
    ## be valid).
    trim1 <- protein_start0 - sum(head(cds_widths, idx1 - 1L))
    trim2 <- cds_cumwidths[[idx2]] - protein_end0 - 1L
    ans <- .trim_first_and_last_ranges(ans, trim1, trim2)

    ## Adjust metadata columns 'protein_start' and 'protein_end' to account
    ## for trimming.
    ans_mcols <- mcols(ans)
    ans_mcols[1L, "protein_start"] <- protein_start
    ans_mcols[nrow(ans_mcols), "protein_end"] <- protein_end
    mcols(ans) <- ans_mcols

    ans
}

### Returns a named GRangesList object parallel to 'x' (names on 'x' are
### propagated).
setMethod("proteinToGenome", "GRangesList",
    function(x, txdb)
    {
        if (!is(x, "IRanges"))
            stop(wmsg("'x' must be an IRanges object or derivative"))
        x_names <- names(x)
        if (is.null(x_names))
            stop(wmsg("'x' must have names"))
        coding_tx_names <- names(txdb)
        if (is.null(coding_tx_names))
            stop(wmsg("'txdb' must have names when it's a GRangesList object"))
        non_coding_idx <- which(!(x_names %in% coding_tx_names))
        if (length(non_coding_idx) != 0L) {
            msg <- .make_bad_names_msg(x_names, non_coding_idx,
                                       what="non-coding transcript name")
            warning(wmsg(msg), immediate.=TRUE)
        }
        x_start <- start(x)
        x_end <- end(x)
        ## TODO: Replace this inefficient lapply-based implementation with
        ## something better.
        ans <- lapply(setNames(seq_along(x), x_names),
            function(i) {
                if (i %in% non_coding_idx)
                    return(GRanges())
                tx_name <- x_names[[i]]
                cds_parts <- txdb[[tx_name]]
                .map_protein_to_cds_parts(x_start[[i]], x_end[[i]], cds_parts)
            }
        )
        GRangesList(ans)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Default proteinToGenome() method
###

### 'txdb' must be a TxDb object or any object that supports transcripts()
### (e.g. EnsDb object).
.check_supplied_tx_names <- function(supplied_tx_names, txdb)
{
    if (is.null(supplied_tx_names))
        stop(wmsg("'x' must have names and they must be transcript names"))
    stopifnot(is.character(supplied_tx_names))
    if (any(supplied_tx_names %in% c(NA_character_, "")))
        stop(wmsg("the names on 'x' cannot contain NAs or empty strings"))
    tx <- transcripts(txdb, columns="tx_name")
    tx_names <- mcols(tx)$tx_name
    bad <- !(supplied_tx_names %in% tx_names)
    if (all(bad))
        stop(wmsg("The names on 'x' must be transcript names present in ",
                  "the supplied ", class(txdb), " object. Note that the ",
                  "transcript names in this object can be obtained/seen ",
                  "with:"),
             "\n    tx <- transcripts(txdb, columns=\"tx_name\")",
             "\n    mcols(tx)$tx_name")
    bad_idx <- which(bad)
    if (length(bad_idx) != 0L) {
        msg <- .make_bad_names_msg(supplied_tx_names, bad_idx,
                                   what="invalid transcript name")
        stop(wmsg(msg))
    }
}

### 'txdb' must be a TxDb object or any object that supports cdsBy()
### (e.g. EnsDb object).
.extract_cds_by_tx <- function(txdb, tx_names)
{
    stopifnot(is.character(tx_names))
    if (!is(txdb, "EnsDb"))
        return(cdsBy(txdb, by="tx", use.names=TRUE))
    ## Should never happen in practice because if 'txdb' is an EnsDb object
    ## then the ensembldb package should be loaded already, and ensembldb
    ## depends on AnnotationFilter.
    if (!requireNamespace("AnnotationFilter", quietly=TRUE))
        stop(wmsg("Couldn't load the AnnotationFilter package. ",
                  "The AnnotationFilter package is needed when ",
                  "calling proteinToGenome() on an EnsDb object. ",
                  "Please install it."))
    filter <- AnnotationFilter::TxIdFilter(tx_names)
    cdsBy(txdb, by="tx", filter=filter)
}

### 'txdb' must be a TxDb object or any object that supports transcripts()
### and cdsBy() (e.g. EnsDb object).
### Returns a named GRangesList object parallel to 'x' (names on 'x' are
### propagated).
.default_proteinToGenome <- function(x, txdb)
{
    if (!is(x, "IRanges"))
        stop(wmsg("'x' must be an IRanges object or derivative"))
    x_names <- names(x)
    .check_supplied_tx_names(x_names, txdb)
    cds_by_tx <- .extract_cds_by_tx(txdb, x_names)
    proteinToGenome(x, cds_by_tx)
}

setMethod("proteinToGenome", "ANY", .default_proteinToGenome)

