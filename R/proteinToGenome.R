### =========================================================================
### proteinToGenome()
### -------------------------------------------------------------------------


### Dispatch is on 2nd argument!
setGeneric("proteinToGenome", signature="txdb",
    function(x, txdb, ...) standardGeneric("proteinToGenome")
)


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
        x_start <- start(x)
        x_end <- end(x)
        ## TODO: Replace this inefficient lapply-based implementation with
        ## something better.
        ans <- lapply(setNames(seq_along(x), x_names),
            function(i) {
                tx_name <- x_names[[i]]
                if (!(tx_name %in% coding_tx_names))
                    return(GRanges())
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
        if (length(bad_idx) == 1L) {
            msg <- "The names on 'x' contain invalid transcript name"
        } else {
            msg <- "The names on 'x' contain invalid transcript names"
            if (length(bad_idx) > 5L) {
                msg <- c(msg, " (showing the first 5 only)")
                bad_idx <- head(bad_idx, n=5L)
            }
        }
        bad_tx_names <- supplied_tx_names[bad_idx]
        stop(wmsg(msg, ": ", paste(bad_tx_names, collapse=", ")))
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

