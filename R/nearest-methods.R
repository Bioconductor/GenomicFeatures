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
        rng <- switch(type,
                      gene=.extractByGeneID(y, id),
                      tx=transcripts(y, "tx_id", filter=list(tx_id=id)),
                      exon=exons(y, "exon_id", filter=list(exon_id=id)),
                      cds=cds(y, "cds_id", filter=list(cds_id=id)))
        if (type != "gene") {
            rng <- .subsetByID(rng, id)
            if (!identical(length(x), length(rng)))
            stop(paste0(type, " regions in annotation cannot be collapsed ",
                        "into a single range"))
        }
        distance(x, rng, ignore.strand=ignore.strand)
    }
)

.extractByGeneID <- function(y, id)
{
    tx <- transcripts(y, "gene_id", filter=list(gene_id=id))
    f <- factor(unlist(tx$gene_id, use.names=FALSE))
    missing <- !id %in% levels(f) 
    if (any(missing))
          stop("'", paste(id[missing], sep=","), "'", " not found in 'y'")

    rng <- unlist(range(split(tx, f)), use.names=FALSE)
    rng <- rng[match(id, levels(f))]
    if (length(rng) != length(id))
        stop("gene regions in annotation 'y' cannot be collapsed ",
             "into a single range")
    rng 
}

.subsetByID <- function(rng, id)
{
    f <- factor(mcols(rng)[,])
    missing <- !id %in% levels(f) 
    if (any(missing))
          stop("'", paste(id[missing], sep=","), "'", " not found in 'y'")
    rng[match(id, levels(f))]
}

