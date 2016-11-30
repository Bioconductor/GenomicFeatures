### =========================================================================
### nearest (and related) methods
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### distance 
###

setMethod("distance", c("GenomicRanges", "TxDb"),
    function(x, y, ignore.strand=FALSE, ..., id, 
             type=c("gene", "tx", "exon", "cds"))
    {
        if (!identical(length(x), length(id)))
            stop("length(id) must equal length(x)")
        if (!is.character(id))
            stop("'id' must be a character")

        if (type == "gene") {
            .extractByGeneID(x, y, ignore.strand, id)
        } else {
            rng <- switch(type,
                          tx=transcripts(y, "tx_id", filter=list(tx_id=id)),
                          exon=exons(y, "exon_id", filter=list(exon_id=id)),
                          cds=cds(y, "cds_id", filter=list(cds_id=id)))
            f <- factor(mcols(rng)[,])
            missing <- !id %in% levels(f) 
            if (any(missing))
                  warning(paste0("id(s): '", paste(unique(id[missing]), 
                                 sep=","), "' were not found in 'y'"))
            ## rep out ranges according to 'id'
            rng <- rng[match(id[!missing], levels(f))]
            ans <- rep(NA_integer_, length(x))
            ans[!missing] <- distance(x[!missing], rng, 
                                      ignore.strand=ignore.strand)
            stopifnot(length(ans) == length(id))
            ans
        }
    }
)

.extractByGeneID <- function(x, y, ignore.strand, id)
{
    tx <- transcriptsBy(y, "gene")
    missing <- !id %in% names(tx) 
    if (any(missing))
          warning(paste0("id(s): '", paste(unique(id[missing]), sep=","), 
                  "' were not found in 'y'"))

    group <- range(tx[names(tx) %in% id], ignore.strand=ignore.strand)
    multiRange <- lengths(group) > 1L
    if (any(multiRange)) {
        warning(paste0("id(s): '", paste(unique(names(multiRange)[multiRange]), 
                       sep=','), 
                       "' could not be collapsed to a single gene region"))
        group <- group[!multiRange]
    }

    valid <- (!id %in% names(multiRange)[multiRange]) & !missing 
    ## rep out ranges according to 'id'
    rng <- unlist(group, use.names=FALSE)
    rng <- rng[match(id[valid], names(group))]
    ans <- rep(NA_integer_, length(x))
    ans[valid] <- distance(x[valid], rng, ignore.strand=ignore.strand)
    stopifnot(length(ans) == length(id))
    ans
}
