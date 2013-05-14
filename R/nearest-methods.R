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
                      gene=.extractGenes(y, id, type),
                      tx=transcripts(y, list(tx_id=id), "tx_id"),
                      exon=exons(y, list(exon_id=id), "exon_id"),
                      cds=cds(y, list(cds_id=id), "cds_id"))
        if (type != "gene")
            .checkID(mcols(rng)[[1]], id, type) 
        if (!identical(length(x), length(rng)))
            stop(paste0(type, " regions in annotation 'y' cannot be collapsed ",
                        "into a single range"))
        distance(x, rng, ignore.strand=ignore.strand)
    }
)

.extractGenes <- function(y, id, type)
{
    tx <- transcripts(y, list(gene_id=id), "gene_id")
    nms <- unlist(tx$gene_id, use.names=FALSE)
    .checkID(nms, id, type)
    skel <- CharacterList(factor(unlist(tx$gene_id, use.names=FALSE)))
    unlist(range(relist(tx, skel)), use.names=FALSE)
}

.checkID <- function(nms, id, type)
{
    missing <- !id %in% nms 
    if (any(missing))
          stop("'", paste(id[missing], sep=","), "'", " not found in 'y'")
}














