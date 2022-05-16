# disjointExons() ---------------------------------------------------------

#' Extract non-overlapping exon parts from an object
#'
#' \code{disjointExons} extracts the non-overlapping exon parts from a
#' \link{TxDb} object or any other supported object.
#'
#' WARNING: \code{disjointExons} is defunct in BioC 3.15. Please use
#' \code{\link{exonicParts}} instead.
#'
#' \code{disjointExons} creates a \link[GenomicRanges]{GRanges} of
#' non-overlapping exon parts with metadata columns of gene_id and exonic_part.
#' Exon parts that overlap more than 1 gene can be dropped with
#' \code{aggregateGenes=FALSE}. When \code{includeTranscripts=TRUE} a
#' \code{tx_name} metadata column is included that lists all transcript names
#' that overlap the exon fragment. This function replaces
#' \code{prepareAnnotationForDEXSeq} in the \pkg{DEXSeq} package.
#'
#' @aliases disjointExons disjointExons,TxDb-method
#' @param x A \link{TxDb} object or any other supported object.
#' @param ...  Arguments to be passed to methods.
#' @param aggregateGenes For \code{disjointExons} : A \code{logical}. When
#' \code{FALSE} (default) exon fragments that overlap multiple genes are
#' dropped.  When \code{TRUE}, all fragments are kept and the \code{gene_id}
#' metadata column includes all gene ids that overlap the exon fragment.
#' @param includeTranscripts For \code{disjointExons} : A \code{logical}. When
#' \code{TRUE} (default) a \code{tx_name} metadata column is included that
#' lists all transcript names that overlap the exon fragment.
#' @return A \link[GenomicRanges]{GRanges} object.
#' @author \code{disjointExons} was originally implemented by Mike Love and
#' Alejandro Reyes and then moved (and adapted) to \pkg{GenomicFeatures} by
#' Valerie Obenchain.
#' @seealso \code{\link{exonicParts}} for an improved version of
#' \code{disjointExons}.
#' @keywords methods
#'
#' @export
setGeneric("disjointExons",
    function(x, ...)
             standardGeneric("disjointExons")
)

setMethod("disjointExons", "TxDb",
    function(x, aggregateGenes=FALSE, includeTranscripts=TRUE, ...)
    {
        msg <- "  disjointExons() is defunct in BioC 3.15. Please use exonicParts() instead."
        .Defunct(msg=wmsg(msg))
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

