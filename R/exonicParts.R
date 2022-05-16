### =========================================================================
### Extraction of exonic and intronic parts
### -------------------------------------------------------------------------
###
### For all functions in this file, 'txdb' must be a TxDb object or any
### object that supports transcripts() and exonsBy() (e.g. EnsDb object).
###


### Works on whatever 'x' can be used as a splitting factor in splitAsList().
### TODO: Rename and move to a more appropriate place (IRanges?)
.rank_in_group <- function(x)
{
    groups <- splitAsList(seq_along(x), x)
    i <- unlist(groups, use.names=FALSE)
    ans <- sequence(lengths(groups))
    ans[i] <- ans
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 3 helper functions used internally by exonicParts() and intronicParts()
###

### Return a GRanges object with 1 range per transcript and metadata columns
### tx_id, tx_name, and gene_id.
### If 'drop.geneless' is FALSE (the default) then the transcripts are
### returned in the same order as with transcripts(), which is expected
### to be by transcript id (tx_id). Otherwise they are ordered first by
### gene id (gene_id), then by transcript id.
tidyTranscripts <- function(txdb, drop.geneless=FALSE)
{
    tx <- transcripts(txdb, columns=c("tx_id", "tx_name", "gene_id"))
    mcols(tx)$gene_id <- as.character(mcols(tx)$gene_id)
    if (drop.geneless) {
        gene_id <- mcols(tx)$gene_id
        tx_id <- mcols(tx)$tx_id
        tx <- tx[order(gene_id, tx_id, na.last=NA)]
    }
    tx
}

### Return a GRangesList object parallel to 'tx_ids'. The supplied 'tx_ids'
### must be a subset of 'mcols(transcripts(txdb))$tx_id'.
.exons_by_txids <- function(txdb, tx_ids)
{
    if (anyDuplicated(tx_ids))
        stop(wmsg("\"transcripts\" method for ", class(txdb), " objects ",
                  "seems broken, sorry"))
    ans <- exonsBy(txdb, by="tx")
    tx_ids <- as.character(tx_ids)
    ans_names <- names(ans)
    if (!identical(tx_ids, ans_names)) {
        m <- match(tx_ids, ans_names)
        if (anyNA(m))
            stop(wmsg("\"exonsBy\" method for ", class(txdb), " objects ",
                      "seems broken, sorry"))
        ans <- ans[m]
    }
    ans
}

### Return a GRanges object with 1 range per exon and metadata columns
### tx_id, tx_name, gene_id, exon_id, exon_name, and exon_rank.
### If 'drop.geneless' is FALSE (the default) then the exons are ordered first
### by transcript id (tx_id), then by exon rank (exon_rank). Otherwise they
### are ordered first by gene id (gene_id), then by transcript id, and then
### by exon rank.
tidyExons <- function(txdb, drop.geneless=FALSE)
{
    tx <- tidyTranscripts(txdb, drop.geneless=drop.geneless)
    ex_by_tx <- .exons_by_txids(txdb, mcols(tx)$tx_id)

    ans <- unlist(ex_by_tx, use.names=FALSE)
    idx <- rep(seq_along(tx), lengths(ex_by_tx))
    mcols(ans) <- cbind(mcols(tx)[idx, , drop=FALSE], mcols(ans))
    ans
}

### Return a GRanges object with 1 range per intron and metadata columns
### tx_id, tx_name, and gene_id.
### If 'drop.geneless' is FALSE (the default) then the introns are ordered
### by transcript id (tx_id). Otherwise they are ordered first by gene id
### (gene_id), then by transcript id.
tidyIntrons <- function(txdb, drop.geneless=FALSE)
{
    tx <- tidyTranscripts(txdb, drop.geneless=drop.geneless)
    ex_by_tx <- .exons_by_txids(txdb, mcols(tx)$tx_id)

    introns_by_tx <- psetdiff(tx, ex_by_tx)

    ans <- unlist(introns_by_tx, use.names=FALSE)
    idx <- rep(seq_along(tx), lengths(introns_by_tx))
    mcols(ans) <- mcols(tx)[idx, , drop=FALSE]
    ans
}

.break_in_parts <- function(x, linked.to.single.gene.only=FALSE,
                               extra_mcol="exonic_part")
{
    ans <- disjoin(x, with.revmap=TRUE)
    revmap <- mcols(ans)$revmap
    ans_mcols <- lapply(mcols(x),
                        function(col) {
                            col <- unique(extractList(col, revmap))
                            col[!is.na(col)]
                        }
                 )
    mcols(ans) <- DataFrame(ans_mcols)
    if (linked.to.single.gene.only) {
        keep_idx <- which(elementNROWS(mcols(ans)$gene_id) == 1L)
        ans <- ans[keep_idx]
        gene_id <- as.character(mcols(ans)$gene_id)
        mcols(ans)$gene_id <- gene_id
        ## Add "exonic_part" or "intronic_part" metadata column for
        ## compatibility with old disjointExons().
        mcols(ans)[[extra_mcol]] <- .rank_in_group(gene_id)
    }
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### exonicParts() and intronicParts()
###

### Return a disjoint and strictly sorted GRanges object with 1 range per
### exonic part and with metadata columns tx_id, tx_name, gene_id, exon_id,
### exon_name, and exon_rank.


#' Extract non-overlapping exonic or intronic parts from a TxDb-like object
#' 
#' \code{exonicParts} and \code{intronicParts} extract the non-overlapping
#' (a.k.a. disjoint) exonic or intronic parts from a \link{TxDb}-like object.
#' 
#' 
#' @aliases tidyTranscripts tidyExons tidyIntrons exonicParts intronicParts
#' @param txdb A \link{TxDb} object, or any \link{TxDb}-like object that
#' supports the \code{\link{transcripts}()} and \code{\link{exonsBy}()}
#' extractors (e.g. an \link[ensembldb]{EnsDb} object).
#' @param linked.to.single.gene.only \code{TRUE} or \code{FALSE}.
#' 
#' If \code{FALSE} (the default), then the disjoint parts are obtained by
#' calling \code{\link[IRanges]{disjoin}()} on all the exons (or introns) in
#' \code{txdb}, including on exons (or introns) not linked to a gene or linked
#' to more than one gene.
#' 
#' If \code{TRUE}, then the disjoint parts are obtained in 2 steps: \enumerate{
#' \item call \code{\link[IRanges]{disjoin}()} on the exons (or introns) linked
#' to \emph{at least one gene},
#' 
#' \item then drop the parts linked to more than one gene from the set of
#' exonic (or intronic) parts obtained previously.  }
#' @param drop.geneless If \code{FALSE} (the default), then all the transcripts
#' (or exons, or introns) get extracted from the \link{TxDb} object.
#' 
#' If \code{TRUE}, then only the transcripts (or exons, or introns) that are
#' linked to a gene get extracted from the \link{TxDb} object.
#' 
#' Note that \code{drop.geneless} also impacts the order in which the features
#' are returned: \itemize{ \item Transcripts: If \code{drop.geneless} is
#' \code{FALSE} then transcripts are returned in the same order as with
#' \code{\link{transcripts}}, which is expected to be by internal transcript id
#' (\code{tx_id}).  Otherwise they are ordered first by gene id
#' (\code{gene_id}), then by internal transcript id.  \item Exons: If
#' \code{drop.geneless} is \code{FALSE} then exons are ordered first by
#' internal transcript id (\code{tx_id}), then by exon rank (\code{exon_rank}).
#' Otherwise they are ordered first by gene id (\code{gene_id}), then by
#' internal transcript id, and then by exon rank.  \item Introns: If
#' \code{drop.geneless} is \code{FALSE} then introns are ordered by internal
#' transcript id (\code{tx_id}).  Otherwise they are ordered first by gene id
#' (\code{gene_id}), then by internal transcript id.  }
#' @return \code{exonicParts} returns a disjoint and strictly sorted
#' \link[GenomicRanges]{GRanges} object with 1 range per exonic part and with
#' metadata columns \code{tx_id}, \code{tx_name}, \code{gene_id},
#' \code{exon_id}, \code{exon_name}, and \code{exon_rank}.  If
#' \code{linked.to.single.gene.only} was set to \code{TRUE}, an additional
#' \code{exonic_part} metadata column is added that indicates the rank of each
#' exonic part within all the exonic parts linked to the same gene.
#' 
#' \code{intronicParts} returns a disjoint and strictly sorted
#' \link[GenomicRanges]{GRanges} object with 1 range per intronic part and with
#' metadata columns \code{tx_id}, \code{tx_name}, and \code{gene_id}.  If
#' \code{linked.to.single.gene.only} was set to \code{TRUE}, an additional
#' \code{intronic_part} metadata column is added that indicates the rank of
#' each intronic part within all the intronic parts linked to the same gene.
#' 
#' \code{tidyTranscripts} returns a \link[GenomicRanges]{GRanges} object with 1
#' range per transcript and with metadata columns \code{tx_id}, \code{tx_name},
#' and \code{gene_id}.
#' 
#' \code{tidyExons} returns a \link[GenomicRanges]{GRanges} object with 1 range
#' per exon and with metadata columns \code{tx_id}, \code{tx_name},
#' \code{gene_id}, \code{exon_id}, \code{exon_name}, and \code{exon_rank}.
#' 
#' \code{tidyIntrons} returns a \link[GenomicRanges]{GRanges} object with 1
#' range per intron and with metadata columns \code{tx_id}, \code{tx_name}, and
#' \code{gene_id}.
#' @note \code{exonicParts} is a replacement for \code{\link{disjointExons}}
#' with the following differences/improvements: \itemize{ \item Argument
#' \code{linked.to.single.gene.only} in \code{exonicParts} replaces argument
#' \code{aggregateGenes} in \code{disjointExons}, but has opposite meaning i.e.
#' \code{exonicParts(txdb, linked.to.single.gene.only=TRUE)} returns the same
#' exonic parts as \code{disjointExons(txdb, aggregateGenes=FALSE)}.
#' 
#' \item Unlike \code{disjointExons(txdb, aggregateGenes=TRUE)},
#' \code{exonicParts(txdb, linked.to.single.gene.only=FALSE)} does NOT discard
#' exon parts that are not linked to a gene.
#' 
#' \item \code{exonicParts} is almost 2x more efficient than
#' \code{disjointExons}.
#' 
#' \item \code{exonicParts} works out-of-the-box on any \link{TxDb}-like object
#' that supports the \code{\link{transcripts}()} and \code{\link{exonsBy}()}
#' extractors (e.g. on an \link[ensembldb]{EnsDb} object).  }
#' @author Hervé Pagès
#' @seealso \itemize{ \item \code{\link[IRanges]{disjoin}} in the \pkg{IRanges}
#' package.
#' 
#' \item \code{\link{transcripts}}, \code{\link{transcriptsBy}}, and
#' \code{\link{transcriptsByOverlaps}}, for extracting genomic feature
#' locations from a \link{TxDb}-like object.
#' 
#' \item \code{\link{transcriptLengths}} for extracting the transcript lengths
#' (and other metrics) from a \link{TxDb} object.
#' 
#' \item \code{\link{extendExonsIntoIntrons}} for extending exons into their
#' adjacent introns.
#' 
#' \item \code{\link{extractTranscriptSeqs}} for extracting transcript (or CDS)
#' sequences from chromosome sequences.
#' 
#' \item \code{\link{coverageByTranscript}} for computing coverage by
#' transcript (or CDS) of a set of ranges.
#' 
#' \item The \link{TxDb} class.  }
#' @keywords manip
#' @examples
#' 
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' 
#' ## ---------------------------------------------------------------------
#' ## exonicParts()
#' ## ---------------------------------------------------------------------
#' 
#' exonic_parts1 <- exonicParts(txdb)
#' exonic_parts1
#' 
#' ## Mapping from exonic parts to genes is many-to-many:
#' gene_id1 <- mcols(exonic_parts1)$gene_id
#' gene_id1  # CharacterList object
#' table(lengths(gene_id1))
#' ## The number of known genes a Human exonic part can be linked to
#' ## varies from 0 to 22!
#' 
#' exonic_parts2 <- exonicParts(txdb, linked.to.single.gene.only=TRUE)
#' exonic_parts2
#' 
#' ## Mapping from exonic parts to genes now is many-to-one:
#' gene_id2 <- mcols(exonic_parts2)$gene_id
#' gene_id2[1:20]  # character vector
#' 
#' ## Select exonic parts for a given gene:
#' exonic_parts2[gene_id2 %in% "643837"]
#' 
#' ## Sanity checks:
#' stopifnot(isDisjoint(exonic_parts1), isStrictlySorted(exonic_parts1))
#' stopifnot(isDisjoint(exonic_parts2), isStrictlySorted(exonic_parts2))
#' stopifnot(all(exonic_parts2 %within% reduce(exonic_parts1)))
#' stopifnot(identical(
#'     lengths(gene_id1) == 1L,
#'     exonic_parts1 %within% exonic_parts2
#' ))
#' 
#' ## ---------------------------------------------------------------------
#' ## intronicParts()
#' ## ---------------------------------------------------------------------
#' 
#' intronic_parts1 <- intronicParts(txdb)
#' intronic_parts1
#' 
#' ## Mapping from intronic parts to genes is many-to-many:
#' mcols(intronic_parts1)$gene_id
#' table(lengths(mcols(intronic_parts1)$gene_id))
#' ## A Human intronic part can be linked to 0 to 22 known genes!
#' 
#' intronic_parts2 <- intronicParts(txdb, linked.to.single.gene.only=TRUE)
#' intronic_parts2
#' 
#' ## Mapping from intronic parts to genes now is many-to-one:
#' class(mcols(intronic_parts2)$gene_id)  # character vector
#' 
#' ## Sanity checks:
#' stopifnot(isDisjoint(intronic_parts1), isStrictlySorted(intronic_parts1))
#' stopifnot(isDisjoint(intronic_parts2), isStrictlySorted(intronic_parts2))
#' stopifnot(all(intronic_parts2 %within% reduce(intronic_parts1)))
#' stopifnot(identical(
#'     lengths(mcols(intronic_parts1)$gene_id) == 1L,
#'     intronic_parts1 %within% intronic_parts2
#' ))
#' 
#' ## ---------------------------------------------------------------------
#' ## Helper functions
#' ## ---------------------------------------------------------------------
#' 
#' tidyTranscripts(txdb)                      # Ordered by 'tx_id'.
#' tidyTranscripts(txdb, drop.geneless=TRUE)  # Ordered first by 'gene_id',
#'                                            # then by 'tx_id'.
#' 
#' tidyExons(txdb)                            # Ordered first by 'tx_id',
#'                                            # then by 'exon_rank'.
#' tidyExons(txdb, drop.geneless=TRUE)        # Ordered first by 'gene_id',
#'                                            # then by 'tx_id',
#'                                            # then by 'exon_rank'.
#' 
#' tidyIntrons(txdb)                          # Ordered by 'tx_id'.
#' tidyIntrons(txdb, drop.geneless=TRUE)      # Ordered first by 'gene_id',
#'                                            # then by 'tx_id'.
#' 
#' @export exonicParts
exonicParts <- function(txdb, linked.to.single.gene.only=FALSE)
{
    if (!isTRUEorFALSE(linked.to.single.gene.only))
        stop("'linked.to.single.gene.only' must be TRUE or FALSE")
    ex <- tidyExons(txdb, drop.geneless=linked.to.single.gene.only)
    .break_in_parts(ex, linked.to.single.gene.only,
                        extra_mcol="exonic_part")
}

### Return a disjoint and strictly sorted GRanges object with 1 range per
### intronic part and with metadata columns tx_id, tx_name, and gene_id.
intronicParts <- function(txdb, linked.to.single.gene.only=FALSE)
{
    if (!isTRUEorFALSE(linked.to.single.gene.only))
        stop("'linked.to.single.gene.only' must be TRUE or FALSE")
    introns <- tidyIntrons(txdb, drop.geneless=linked.to.single.gene.only)
    .break_in_parts(introns, linked.to.single.gene.only,
                             extra_mcol="intronic_part")
}

