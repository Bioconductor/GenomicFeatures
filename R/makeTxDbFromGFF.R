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
                                taxonomyId=NA, miRBaseBuild=NA)
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
    if(is.na(taxonomyId)){
        taxonomyId <- GenomeInfoDb:::.taxonomyId(organism)
    }else{
        GenomeInfoDb:::.checkForAValidTaxonomyId(taxonomyId)
    }
    metadata <- data.frame(
                   name=c("Data source",
                          "Organism",
                          "Taxonomy ID",
                          "miRBase build ID"),
                   value=c(dataSource, organism,
                     taxonomyId,
                     miRBaseBuild)
                   )
    message("metadata: OK")
    metadata
}

.detect_file_format <- function(file)
{
    if (is(file, "connection"))
        file <- summary(file)$description
    if (isSingleString(file)) {
        file2 <- try(FileForFormat(file), silent=TRUE)
        if (inherits(file2, "try-error"))
            return(tools::file_ext(file))
        file <- file2
    }
    if (is(file, "RTLFile")) {
        if (is(file, "GFF3File"))
            return("gff3")
        if (is(file, "GTFFile"))
            return("gtf")
        desc <- rtracklayer:::resourceDescription(file)
        if (is(file, "CompressedFile")) {
            ## Mimic what import,CompressedFile,missing,missing method does.
            desc <- tools::file_path_sans_ext(desc)
        }
        format <- tools::file_ext(desc)
        return(format)
    }
    stop(wmsg("Invalid 'file'. Must be a path to a file, or an URL, ",
              "or a connection object, or a GFF3File or GTFFile object."))
}

### Based on makeTxDbFromGRanges().
makeTxDbFromGFF <- function(file,
                            format=c("auto", "gff3", "gtf"),
                            dataSource=NA,
                            organism=NA,
                            taxonomyId=NA,
                            circ_seqs=DEFAULT_CIRC_SEQS,
                            chrominfo=NA,
                            miRBaseBuild=NA,
                            exonRankAttributeName=NA,     # defunct
                            gffGeneIdAttributeName=NA,    # defunct
                            useGenesAsTranscripts=FALSE,  # defunct
                            gffTxName="mRNA",             # defunct
                            species=NA)                   # defunct
{
    ## Raise error for a bunch of args that are defunct.
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
    if (format == "auto") {
        format <- .detect_file_format(file)
        if (!(format %in% c("gff3", "gff", "gtf")))
            stop(wmsg("Cannot detect whether 'file' is a GFF3 or GTF file. ",
                      "Please use the 'format' argument to specify the ",
                      "format (\"gff3\" or \"gtf\")."))
    }
    gr <- import(file, format=format, feature.type=GFF_FEATURE_TYPES)
    gr <- .set_seqinfo(gr, chrominfo, circ_seqs)
    metadata <- .prepareGFFMetadata(file, dataSource, organism, taxonomyId,
                                    miRBaseBuild)
    makeTxDbFromGRanges(gr, metadata=metadata)
}

makeTranscriptDbFromGFF <- function(...) .Defunct("makeTxDbFromGFF")

