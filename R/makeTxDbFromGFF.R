### =========================================================================
### makeTxDbFromGFF()
### -------------------------------------------------------------------------


.set_seqinfo <- function(gr, chrominfo=NA, circ_seqs=DEFAULT_CIRC_SEQS)
{
    if (identical(chrominfo, NA)) {
        seqlevels <- seqlevels(gr)
        seqlevels[rankSeqlevels(seqlevels)] <- seqlevels
        seqlevels(gr) <- seqlevels
        isCircular(gr) <- matchCircularity(seqlevels(gr), circ_seqs)
        return(gr)
    }
    if (!identical(circ_seqs, DEFAULT_CIRC_SEQS))
        stop(wmsg("only one of 'chrominfo' and 'circ_seqs' can be specified"))
    if (!is.data.frame(chrominfo))
        stop(wmsg("'chrominfo' must be a data.frame or NA"))
    if (!identical(colnames(chrominfo),
                   c("chrom", "length", "is_circular")))
        stop(wmsg("'chrominfo' must have three  columns that correpspond ",
                  "to 'chrom', 'length', and 'is_circular' and are named ",
                  "accordingly"))
    seqlevels(gr) <- as.character(chrominfo$chrom)
    seqlengths(gr) <- chrominfo$length
    isCircular(gr) <- chrominfo$is_circular
    gr
}

.prepareGFFMetadata <- function(file, dataSource=NA, organism=NA,
                                      miRBaseBuild=NA)
{
    message("Prepare the 'metadata' data frame ... ",
            appendLF=FALSE)
    if (!isSingleStringOrNA(dataSource))
        stop("'dataSource' must be a a single string or NA")
    if (!isSingleStringOrNA(organism))
        stop("'organism' must be a a single string or NA")
    if (!isSingleStringOrNA(miRBaseBuild))
        stop("'miRBaseBuild' must be a a single string or NA")
    if (identical(dataSource, NA)) {
        if (is.character(file)) {
            dataSource <- file
        } else {
            dataSource <- showConnections(all=TRUE)[as.character(file),
                                                    "description"]
        }
    }
    metadata <- data.frame(
                   name=c("Data source",
                          "Organism",
                          "miRBase build ID"),
                   value=c(dataSource, organism, miRBaseBuild)
                   )
    message("metadata: OK")
    metadata
}

### Based on makeTxDbFromGRanges().
makeTxDbFromGFF <- function(file,
                            format=c("gff3", "gtf"),
                            exonRankAttributeName=NA,     # deprecated
                            gffGeneIdAttributeName=NA,    # deprecated
                            chrominfo=NA,
                            dataSource=NA,
                            organism=NA,
                            circ_seqs=DEFAULT_CIRC_SEQS,
                            miRBaseBuild=NA,
                            useGenesAsTranscripts=FALSE,  # deprecated
                            gffTxName="mRNA",             # deprecated
                            species=NA)                   # deprecated
{
    if (!identical(exonRankAttributeName, NA))
        .Deprecated(msg="'exonRankAttributeName' is ignored and deprecated")
    if (!identical(gffGeneIdAttributeName, NA))
        .Deprecated(msg="'gffGeneIdAttributeName' is ignored and deprecated")
    if (!identical(useGenesAsTranscripts, FALSE))
        .Deprecated(msg="'useGenesAsTranscripts' is ignored and deprecated")
    if (!identical(gffTxName, "mRNA"))
        .Deprecated(msg="'gffTxName' is ignored and deprecated")
    if (!identical(species, NA)) {
        if (!identical(organism, NA))
            stop("only one of 'organism' or 'species' can be specified, ",
                 "but not both")
        msg <- c("The 'species' argument is deprecated. ",
                 "Please use 'organism' instead.")
        .Deprecated(msg=msg)
        organism <- species
    }

    format <- match.arg(format)
    gr <- import(file, format=format)
    gr <- .set_seqinfo(gr, chrominfo, circ_seqs)
    metadata <- .prepareGFFMetadata(file, dataSource, organism, miRBaseBuild)
    makeTxDbFromGRanges(gr, metadata=metadata)
}

makeTranscriptDbFromGFF <- function(...)
{
    .Deprecated("makeTxDbFromGFF")
    makeTxDbFromGFF(...)
}

