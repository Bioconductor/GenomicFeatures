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
                          "TaxID",
                          "miRBase build ID"),
                   value=c(dataSource, organism,
                     GenomeInfoDb:::.taxonomyId(organism),
                     miRBaseBuild)
                   )
    message("metadata: OK")
    metadata
}

### Based on makeTxDbFromGRanges().
makeTxDbFromGFF <- function(file,
                            format=c("gff3", "gtf"),
                            exonRankAttributeName=NA,     # defunct
                            gffGeneIdAttributeName=NA,    # defunct
                            chrominfo=NA,
                            dataSource=NA,
                            organism=NA,
                            circ_seqs=DEFAULT_CIRC_SEQS,
                            miRBaseBuild=NA,
                            useGenesAsTranscripts=FALSE,  # defunct
                            gffTxName="mRNA",             # defunct
                            species=NA)                   # defunct
{
    if (!identical(exonRankAttributeName, NA))
        .Defunct(msg="the 'exonRankAttributeName' argument is defunct")
    if (!identical(gffGeneIdAttributeName, NA))
        .Defunct(msg="the 'gffGeneIdAttributeName' argument is defunct")
    if (!identical(useGenesAsTranscripts, FALSE))
        .Defunct(msg="the 'useGenesAsTranscripts' argument is defunct")
    if (!identical(gffTxName, "mRNA"))
        .Defunct(msg="the 'gffTxName' argument is defunct")
    if (!identical(species, NA))
        .Defunct(msg=wmsg("The 'species' argument is defunct. ",
                          "Please use 'organism' instead."))

    format <- match.arg(format)
    gr <- import(file, format=format, feature.type=GFF_FEATURE_TYPES)
    gr <- .set_seqinfo(gr, chrominfo, circ_seqs)
    metadata <- .prepareGFFMetadata(file, dataSource, organism, miRBaseBuild)
    makeTxDbFromGRanges(gr, metadata=metadata)
}

makeTranscriptDbFromGFF <- function(...) .Defunct("makeTxDbFromGFF")

