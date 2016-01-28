### =========================================================================
### disjointExons()
### -------------------------------------------------------------------------


setGeneric("disjointExons", 
    function(x, ...) 
             standardGeneric("disjointExons")
)

setMethod("disjointExons", "TxDb", 
    function(x, aggregateGenes=FALSE, includeTranscripts=TRUE, ...) 
    {
        exonsByGene <- exonsBy(x, by="gene")
        exonicParts <- disjoin(unlist(exonsByGene, use.names=FALSE))

        if (aggregateGenes) {
            foGG <- findOverlaps(exonsByGene, exonsByGene)
            aggregateNames <- .listNames(names(exonsByGene), as.list(foGG)) 
            foEG <- findOverlaps(exonicParts, exonsByGene, select="first")
            gene_id <- aggregateNames[foEG]
            pasteNames <- .pasteNames(names(exonsByGene), as.list(foGG))[foEG]
            orderByGeneName <- order(pasteNames) 
            exonic_rle <- runLength(Rle(pasteNames[orderByGeneName])) 
        } else {
            ## drop exonic parts that overlap > 1 gene
            foEG <- findOverlaps(exonicParts, exonsByGene)
            idxList <- as.list(foEG)
            if (any(keep <- countQueryHits(foEG) == 1)) {
                idxList <- idxList[keep]
                exonicParts <- exonicParts[keep]
            }
            gene_id <- .listNames(names(exonsByGene), idxList)
            orderByGeneName <- order(unlist(gene_id, use.names=FALSE)) 
            exonic_rle <- runLength(Rle(unlist(gene_id[orderByGeneName],
                                               use.names=FALSE))) 
        }
        values <- DataFrame(gene_id)

        if (includeTranscripts) {
           exonsByTx <- exonsBy(x, by="tx", use.names=TRUE )
           foET <- findOverlaps(exonicParts, exonsByTx)
           values$tx_name <- .listNames(names(exonsByTx), as.list(foET))
        } 
        mcols(exonicParts) <- values 
        exonicParts <- exonicParts[orderByGeneName]
        exonic_part <- unlist(lapply(exonic_rle, seq_len), use.names=FALSE)
        exonicParts$exonic_part <- exonic_part
        exonicParts
    }
)

## returns a character vector the same length as indexList
.pasteNames <- function(names, indexList)
 {
    nm <- names[unlist(indexList, use.names=FALSE)]
    rl <- relist(nm, indexList)
    el <- elementNROWS(indexList) > 1
    rl[el] <- base::lapply(rl[el], base::paste, collapse="+")
    unlist(rl, use.names=FALSE)
}

## returns a CharacterList the same length as indexList
.listNames <- function(names, indexList)
{
    nm <- names[unlist(indexList, use.names=FALSE)]
    unname(CharacterList(relist(nm, indexList)))
}

