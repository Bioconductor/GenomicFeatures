### =========================================================================
### extractUpstreamSeqs()
### -------------------------------------------------------------------------


### Dispatch is on the 2nd argument!
setGeneric("extractUpstreamSeqs", signature="genes",
    function(x, genes, width=1000, ...) standardGeneric("extractUpstreamSeqs")
)

### Will work on any object 'x' for which seqinfo() and getSeq() are defined
### e.g. BSgenome, FaFile, TwoBitFile, etc...
setMethod("extractUpstreamSeqs", "GenomicRanges",
    function(x, genes, width=1000)
    {
        seqinfo(genes) <- merge(seqinfo(genes), seqinfo(x))
        upstream <- trim(suppressWarnings(flank(genes, width=width)))
        ans <- getSeq(x, upstream)

        ## Add metada columns to 'ans'.
        gene_seqnames <- seqnames(genes)
        gene_strand <- strand(genes)
        idx1 <- which(gene_strand != "-")
        idx2 <- which(gene_strand == "-")
        gene_TSS <- integer(length(genes))
        gene_TSS[idx1] <- start(genes)[idx1]
        gene_TSS[idx2] <- end(genes)[idx2]
        ans_mcols <- DataFrame(gene_seqnames=gene_seqnames,
                               gene_strand=gene_strand,
                               gene_TSS=gene_TSS)
        mcols(ans) <- ans_mcols
        ans
    }
)

setMethod("extractUpstreamSeqs", "TxDb",
    function(x, genes, width=1000, exclude.seqlevels=NULL)
    {
        genes <- sort(genes(genes))
        ## 'genes' is now a GRanges object.
        if (!is.null(exclude.seqlevels)) {
            if (!is.character(exclude.seqlevels))
                stop("'exclude.seqlevels' must be NULL or a character vector")
            idx <- match(exclude.seqlevels, seqlevels(genes))
            if (any(is.na(idx)))
                stop("'exclude.seqlevels' contains invalid seqlevels")
            seqlevels(genes, force=TRUE) <- seqlevels(genes)[-idx]
        }
        callGeneric(x, genes, width=width)
    }
)

### 'genes' is assumed to contain transcripts grouped by gene e.g. as returned
### by transcriptsBy(..., by="gene").
setMethod("extractUpstreamSeqs", "GRangesList",
    function(x, genes, width=1000)
    {
        stop("NOT READY YET, SORRY!")
    }
)

