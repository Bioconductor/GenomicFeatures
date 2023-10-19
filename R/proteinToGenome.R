### =========================================================================
### proteinToGenome()
### -------------------------------------------------------------------------


### Dispatch is on 2nd argument!
setGeneric("proteinToGenome", signature="db",
    function(x, db, ...) standardGeneric("proteinToGenome")
)


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

### 'cds' must be a GRanges object representing the CDS parts of a given
### transcript/protein.
### 'protein_start' and 'protein_end' must be protein-relative coordinates
### i.e. coordinates (counted in Amino Acids) relative to the protein
### associated with 'cds'.
### Returns a GRanges object.
.map_protein_to_cds <- function(protein_start, protein_end, cds)
{
    stopifnot(isSingleNumber(protein_start),
              isSingleNumber(protein_end),
              protein_start <= protein_end,
              is(cds, "GRanges"))
    nparts <- length(cds)
    cds_widths <- width(cds)
    stopifnot(nparts >= 1L, all(cds_widths >= 1L))
    protein_start <- as.integer(protein_start)
    protein_end <- as.integer(protein_end)
    cds_cumwidths <- cumsum(cds_widths)

    ## Add metadata columns 'protein_start' and 'protein_end' to 'cds'.
    protein_ranges <- .make_protein_ranges_from_cumwidths(cds_cumwidths)
    protein_ranges <- DataFrame(protein_start=start(protein_ranges),
                                protein_end=end(protein_ranges))
    mcols(cds) <- cbind(mcols(cds), protein_ranges)

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
    ans <- cds[idx1:idx2]

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
    function(x, db)
    {
        if (!is(x, "IRanges"))
            stop(wmsg("'x' must be an IRanges object or derivative"))
        x_names <- names(x)
        if (is.null(x_names))
            stop(wmsg("'x' must have names"))
        coding_tx_names <- names(db)
        if (is.null(coding_tx_names))
            stop(wmsg("'db' must have names when it's a GRangesList object"))
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
                cds <- db[[tx_name]]
                .map_protein_to_cds(x_start[[i]], x_end[[i]], cds)
            }
        )
        GRangesList(ans)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Default proteinToGenome() method
###

### 'db' must be a TxDb object or any object that supports transcripts()
### (e.g. EnsDb object).
.check_supplied_tx_names <- function(supplied_tx_names, db)
{
    if (is.null(supplied_tx_names))
        stop(wmsg("'x' must have names and they must be transcript names"))
    stopifnot(is.character(supplied_tx_names))
    if (any(supplied_tx_names %in% c(NA_character_, "")))
        stop(wmsg("the names on 'x' cannot contain NAs or empty strings"))
    tx <- transcripts(db, columns="tx_name")
    tx_names <- mcols(tx)$tx_name
    bad <- !(supplied_tx_names %in% tx_names)
    if (all(bad))
        stop(wmsg("The names on 'x' must be transcript names present in ",
                  "the supplied ", class(db), " object. Note that the ",
                  "transcript names in this object can be obtained/seen ",
                  "with:"),
             "\n    tx <- transcripts(db, columns=\"tx_name\")",
             "\n    mcols(tx)$tx_name")
    bad_idx <- which(bad)
    if (length(bad_idx) != 0L) {
        msg <- .make_bad_names_msg(supplied_tx_names, bad_idx,
                                   what="invalid transcript name")
        stop(wmsg(msg))
    }
}

### 'db' must be a TxDb object or any object that supports cdsBy()
### (e.g. EnsDb object).
.extract_cds_by_tx <- function(db, tx_names)
{
    stopifnot(is.character(tx_names))
    if (!is(db, "EnsDb"))
        return(cdsBy(db, by="tx", use.names=TRUE))
    ## Should never happen in practice because if 'db' is an EnsDb object
    ## then the ensembldb package should be loaded already, and ensembldb
    ## depends on AnnotationFilter.
    if (!requireNamespace("AnnotationFilter", quietly=TRUE))
        stop(wmsg("Couldn't load the AnnotationFilter package. ",
                  "The AnnotationFilter package is needed when ",
                  "calling proteinToGenome() on an EnsDb object. ",
                  "Please install it."))
    filter <- AnnotationFilter::TxIdFilter(tx_names)
    cdsBy(db, by="tx", filter=filter)
}

### 'db' must be a TxDb object or any object that supports transcripts()
### and cdsBy() (e.g. EnsDb object).
### Returns a named GRangesList object parallel to 'x' (names on 'x' are
### propagated).
.default_proteinToGenome <- function(x, db)
{
    if (!is(x, "IRanges"))
        stop(wmsg("'x' must be an IRanges object or derivative"))
    x_names <- names(x)
    .check_supplied_tx_names(x_names, db)
    cds_by_tx <- .extract_cds_by_tx(db, x_names)
    proteinToGenome(x, cds_by_tx)
}

setMethod("proteinToGenome", "ANY", .default_proteinToGenome)

