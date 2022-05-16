### =========================================================================
### Map internal ids to external names for a given feature type.
### -------------------------------------------------------------------------




#' Map internal ids to external names for a given feature type
#' 
#' Utility function for retrieving the mapping from the internal ids to the
#' external names of a given feature type.
#' 
#' Transcripts, exons and CDS in a \link{TxDb} object are stored in seperate
#' tables where the primary key is an integer called \emph{feature internal
#' id}. This id is stored in the \code{"tx_id"} column for transcripts, in the
#' \code{"exon_id"} column for exons, and in the \code{"cds_id"} column for
#' CDS.  Unlike other commonly used ids like Entrez Gene IDs or Ensembl IDs,
#' this internal id was generated at the time the \link{TxDb} object was
#' created and has no meaning outside the scope of this object.
#' 
#' The \code{id2name} function can be used to translate this internal id into a
#' more informative id or name called \emph{feature external name}. This name
#' is stored in the \code{"tx_name"} column for transcripts, in the
#' \code{"exon_name"} column for exons, and in the \code{"cds_name"} column for
#' CDS.
#' 
#' Note that, unlike the feature internal id, the feature external name is not
#' guaranteed to be unique or even defined (the column can contain \code{NA}s).
#' 
#' @param txdb A \link{TxDb} object.
#' @param feature.type The feature type for which the mapping must be
#' retrieved.
#' @return A named character vector where the names are the internal ids and
#' the values the external names.
#' @author Hervé Pagès
#' @seealso \itemize{ \item \code{\link{transcripts}},
#' \code{\link{transcriptsBy}}, and \code{\link{transcriptsByOverlaps}}, for
#' how to extract genomic features from a \link{TxDb} object.  \item The
#' \link{TxDb} class.  }
#' @examples
#' 
#' txdb1_file <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'                           package="GenomicFeatures")
#' txdb1 <- loadDb(txdb1_file)
#' id2name(txdb1, feature.type="tx")[1:4]
#' id2name(txdb1, feature.type="exon")[1:4]
#' id2name(txdb1, feature.type="cds")[1:4]
#' 
#' txdb2_file <- system.file("extdata", "Biomart_Ensembl_sample.sqlite",
#'                           package="GenomicFeatures")
#' txdb2 <- loadDb(txdb2_file)
#' id2name(txdb2, feature.type="tx")[1:4]
#' id2name(txdb2, feature.type="exon")[1:4]
#' id2name(txdb2, feature.type="cds")[1:4]
#' 
#' @export id2name
id2name <- function(txdb, feature.type=c("tx", "exon", "cds"))
{
    if (!is(txdb, "TxDb"))
        stop("'txdb' must be a TxDb object")
    feature.type <- match.arg(feature.type)
    table <- switch(feature.type, tx="transcript", exon="exon", cds="cds")
    columns <- TXDB_table_columns(table)[c("id", "name")]
    df <- TxDb_SELECT_from_INNER_JOIN(txdb, table, columns)
    ans <- df[[columns[2L]]]
    names(ans) <- as.character(df[[columns[1L]]])
    ans
}

