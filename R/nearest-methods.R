### =========================================================================
### nearest (and related) methods
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### distance 
###

setMethod("distance", c("GenomicRanges", "TranscriptDb"),
    function(x, y, ignore.strand=FALSE, ..., id, 
             type=c("gene", "tx", "exon", "cds"))
    {
        if (!identical(length(x), length(id)))
            stop("length(id) must equal length(x)")
        if (!is.character(id))
            stop("'id' must be a character")
        rng <- switch(type,
                      gene=.extractByGeneID(y, id),
                      tx=transcripts(y, list(tx_id=id), "tx_id"),
                      exon=exons(y, list(exon_id=id), "exon_id"),
                      cds=cds(y, list(cds_id=id), "cds_id"))
        if (type != "gene")
            rng <- .subsetByID(rng, id)
        if (!identical(length(x), length(rng)))
            stop("length(x) != length(subject)")
        distance(x, rng, ignore.strand=ignore.strand)
    }
)

.extractByGeneID <- function(y, id)
{
    tx <- transcripts(y, list(gene_id=id), "gene_id")
    f <- factor(unlist(tx$gene_id, use.names=FALSE))
    missing <- !id %in% levels(f) 
    if (any(missing))
          stop("'", paste(id[missing], sep=","), "'", " not found in 'y'")

    rng <- range(split(tx, f))
    rng[match(levels(f), id)]
    if (any(elementLengths(rng) > 1L))
        stop("gene regions in annotation 'y' cannot be collapsed ",
             "into a single range")
    unlist(rng, use.names=TRUE)
}

.subsetByID <- function(rng, id)
{
    f <- factor(mcols(rng)[,])
    missing <- !id %in% levels(f) 
    if (any(missing))
          stop("'", paste(id[missing], sep=","), "'", " not found in 'y'")
    rng[match(id, levels(f))]
}

