### =========================================================================
### makeTranscriptDbFromUCSC()
### -------------------------------------------------------------------------


### Lookup between UCSC tables and tracks in the "Genes and Gene Prediction"
### group that are compatible with makeTranscriptDbFromUCSC(). A table is
### compatible if it has the following cols: name, chrom, strand, txStart,
### txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds.
### Note that, from a strictly technical point of view, the name and
### exonCount cols are not required (i.e. .makeTranscriptDbFromUCSCTxTable()
### could easily be modified to work even when they are missing).
.SUPPORTED_UCSC_TABLES <- c(
  ## tablename (unique key)    track             subtrack
  "knownGene",                 "UCSC Genes",     NA,
  "knownGeneOld3",             "Old UCSC Genes", NA,
  "wgEncodeGencodeManualRel2", "Gencode Genes",  "Genecode Manual",
  "wgEncodeGencodeAutoRel2",   "Gencode Genes",  "Genecode Auto",
  "wgEncodeGencodePolyaRel2",  "Gencode Genes",  "Genecode PolyA",
  "ccdsGene",                  "Consensus CDS",  NA, 
  "refGene",                   "RefSeq Genes",   NA, 
  "xenoRefGene",               "Other RefSeq",   NA, 
  "vegaGene",                  "Vega Genes",     "Vega Protein Genes", 
  "vegaPseudoGene",            "Vega Genes",     "Vega Pseudogenes", 
  "ensGene",                   "Ensembl Genes",  NA, 
  "acembly",                   "AceView Genes",  NA, 
  "sibGene",                   "SIB Genes",      NA, 
  "nscanPasaGene",             "N-SCAN",         "N-SCAN PASA-EST",
  "nscanGene",                 "N-SCAN",         "N-SCAN", 
  "sgdGene",                   "SGD Genes",      NA, 
  "sgpGene",                   "SGP Genes",      NA,
  "geneid",                    "Geneid Genes",   NA, 
  "genscan",                   "Genscan Genes",  NA, 
  "exoniphy",                  "Exoniphy",       NA, 
  "augustusHints",             "Augustus",       "Augustus Hints", 
  "augustusXRA",               "Augustus",       "Augustus De Novo", 
  "augustusAbinitio",          "Augustus",       "Augustus Ab Initio",
  "acescan",                   "ACEScan",        NA
)

supportedUCSCtables <- function()
{
    mat <- matrix(.SUPPORTED_UCSC_TABLES, ncol=3, byrow=TRUE)
    colnames(mat) <- c("tablename", "track", "subtrack")
    data.frame(track=mat[ , "track"], subtrack=mat[ , "subtrack"],
               row.names=mat[ , "tablename"],
               stringsAsFactors=FALSE)
}

### The table names above (unique key) must be used to name the top-level
### elements of the list below. If no suitable tx_name-to-gene_id mapping is
### available in the UCSC database for a supported table, then there is no
### entry in the list below for this table and makeTranscriptDbFromUCSC()
### will leave the gene table empty.
.UCSC_TXNAME2GENEID_MAPDEFS <- list(
    knownGene=list(
        L2Rchain=list(
            c(tablename="knownToLocusLink",
              Lcolname="name",
              Rcolname="value")
        ),
        gene_id_type="Entrez Gene ID"
    ),
    wgEncodeGencodeManualRel2=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeClassesRel2",
              Lcolname="name",
              Rcolname="geneId")
        ),
        gene_id_type="HAVANA Pseudogene ID"
    ),
    wgEncodeGencodeAutoRel2=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeClassesRel2",
              Lcolname="name",
              Rcolname="geneId")
        ),
        gene_id_type="HAVANA Pseudogene ID"
    ),
    wgEncodeGencodePolyaRel2=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeClassesRel2",
              Lcolname="name",
              Rcolname="geneId")
        ),
        gene_id_type="HAVANA Pseudogene ID"
    ),
    refGene=list(
        L2Rchain=list(
            c(tablename="refLink",
              Lcolname="mrnaAcc",
              Rcolname="locusLinkId")
        ),
        gene_id_type="Entrez Gene ID"
    ),
    vegaGene=list(
        L2Rchain=list(
            c(tablename="vegaGtp",
              Lcolname="transcript",
              Rcolname="gene")
        ),
        gene_id_type="HAVANA Pseudogene ID"
    ),
    vegaPseudoGene=list(
        L2Rchain=list(
            c(tablename="vegaGtp",
              Lcolname="transcript",
              Rcolname="gene")
        ),
        gene_id_type="HAVANA Pseudogene ID"
    ),
    ensGene=list(
        L2Rchain=list(
            c(tablename="ensGtp",
              Lcolname="transcript",
              Rcolname="gene")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    sgdGene=list(
        L2Rchain=list(
            c(tablename="sgdIsoforms",
              Lcolname="transcript",
              Rcolname="clusterId"),
            c(tablename="sgdCanonical",
              Lcolname="clusterId",
              Rcolname="transcript")
        ),
        gene_id_type="ID of canonical transcript in cluster"
    )
)

### Returns a named list with 2 elements:
###   $genes: data.frame with tx_name and gene_id cols;
###   $gene_id_type: single string.
.fetchTxName2GeneIdMappingFromUCSC <- function(session, track, Ltablename)
{
    mapdef <- .UCSC_TXNAME2GENEID_MAPDEFS[[Ltablename]]
    if (is.null(mapdef))
        return(list(genes=NULL, gene_id_type="no gene ids"))
    nlink <- length(mapdef$L2Rchain)
    for (i in seq_len(nlink)) {
        L2Rlink <- mapdef$L2Rchain[[i]]
        tablename <- L2Rlink[["tablename"]]
        Lcolname <- L2Rlink[["Lcolname"]]
        Rcolname <- L2Rlink[["Rcolname"]]
        query <- ucscTableQuery(session, track, table=tablename)
        ucsc_table <- getTable(query)
        if (!all(hasCol(ucsc_table, c(Lcolname, Rcolname))))
            stop("expected cols \"", Lcolname, "\" or/and \"",
                 Rcolname, "\" not found in table ", tablename)
        Lcol <- ucsc_table[[Lcolname]]
        Rcol <- ucsc_table[[Rcolname]]
        if (!is.character(Lcol))
            Lcol <- as.character(Lcol)
        if (!is.character(Rcol))
            Rcol <- as.character(Rcol)
        if (i == 1L) {
            tmp <- data.frame(Lcol=Lcol, Rcol=Rcol)
        } else {
            name2val <- Rcol
            names(name2val) <- Lcol
            tmp <- joinDataFrameWithName2Val(tmp, "Rcol", name2val, "Rcol")
        }
    }
    genes <- data.frame(tx_name=tmp$Lcol, gene_id=tmp$Rcol)
    gene_id_type <- mapdef$gene_id_type
    list(genes=genes, gene_id_type=gene_id_type)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the 'transcripts' data frame from UCSC table.
###

.extractTranscriptsFromUCSCTxTable <- function(ucsc_txtable)
{
    transcripts <- data.frame(
        tx_id=seq_len(nrow(ucsc_txtable)),
        tx_name=ucsc_txtable$name,
        tx_chrom=ucsc_txtable$chrom,
        tx_strand=ucsc_txtable$strand,
        tx_start=ucsc_txtable$txStart + 1L,
        tx_end=ucsc_txtable$txEnd
    )
    transcripts
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the 'splicings' data frame from UCSC table.
###
### Exon starts and ends are multi-valued fields (comma-separated) that
### need to be expanded.

.makeExonRank <- function(exonCount, exonStrand)
{
    ans <- lapply(seq_len(length(exonCount)),
        function(i)
        {
            if (exonStrand[i] == "+")
                seq_len(exonCount[i])
            else
                (exonCount[i]):1L
        }
    )
    unlist(ans)
}

### Returns a named list with 2 elements. Each element is itself a list of
### integer vectors with no NAs. The 2 elements have the same "shape".
.extractExonLocsFromUCSCTxTable <- function(ucsc_txtable)
{
    exon_count <- ucsc_txtable$exonCount
    exon_start <- strsplitAsListOfIntegerVectors(ucsc_txtable$exonStarts)
    if (!identical(elementLengths(exon_start), exon_count))
        stop("UCSC data anomaly: 'ucsc_txtable$exonStarts' ",
             "inconsistent with 'ucsc_txtable$exonCount'")
    exon_end <- strsplitAsListOfIntegerVectors(ucsc_txtable$exonEnds)
    if (!identical(elementLengths(exon_end), exon_count))
        stop("UCSC data anomaly: 'ucsc_txtable$exonEnds' ",
             "inconsistent with 'ucsc_txtable$exonCount'")
    return(list(start=exon_start, end=exon_end))
}

### 'cdsStart0', 'cdsEnd1': single integers (resp. 0-based and 1-based).
### 'exon_start0', 'exon_end1': integer vectors of equal lengths (resp.
### 0-based and 1-based) and with no NAs.
### Returns a list with 2 elements, each of them being an integer vector of
### the same length as 'exon_start0' (or 'exon_end1') that may contain NAs.
### Notes:
###   (1) In refGene table, transcript NM_001146685: cds cumulative length is
###       not a multiple of 3:
###                 name chrom strand txStart   txEnd cdsStart  cdsEnd
###         NM_001146685  chr1      + 1351370 1353029  1351370 1353029
###         exonCount       exonStarts         exonEnds id   name2
###                 2 1351370,1352796, 1351628,1353029,  0 TMEM88B     
###         cdsStartStat cdsEndStat exonFrames
###                 cmpl     incmpl       0,0,
###       --> cds lengths: 1351628 - 1351370 -> 258
###                        1353029 - 1352796 -> 233
###       --> cds cum length: 491
###       Note that the cds end is marked as "incomplete" (see the cdsEndStat
###       col) which, according to UCSC, means that "the CDS is NOT completely
###       contained in the alignment at this end". See this post on the Genome
###       mailing list for more information:
###       https://lists.soe.ucsc.edu/pipermail/genome/2005-December/009184.html
###       Note that the post is about the Gencode Genes. Is it reasonable to
###       assume that this applies to RefSeq Genes too?
###   (2) Same thing in ensGene table, transcript ENST00000371841.
###   (3) All transcripts in knowGene table have a cds cumulative length that
###       is a multiple of 3.
### TODO: Investigate (1) and (2).
.extractUCSCCdsStartEnd <- function(cdsStart0, cdsEnd1,
                                    exon_start0, exon_end1, tx_name)
{
    cds_start0 <- cds_end1 <- integer(length(exon_start0))
    cds_start0[] <- NA_integer_
    cds_end1[] <- NA_integer_
    if (cdsStart0 >= cdsEnd1)
        return(list(cds_start0, cds_end1))
    first_exon_with_cds <- which(exon_start0 <= cdsStart0
                                 & cdsStart0 < exon_end1)
    if (length(first_exon_with_cds) != 1L)
        stop("UCSC data ambiguity in transcript ", tx_name,
             ": cannot determine first exon with cds ('cdsStart' ",
             "falls in 0 or more than 1 exon)")
    last_exon_with_cds <- which(exon_start0 < cdsEnd1
                                & cdsEnd1 <= exon_end1)
    if (length(last_exon_with_cds) != 1L)
        stop("UCSC data ambiguity in transcript ", tx_name,
             ": cannot determine last exon with cds ('cdsEnd' ",
             "falls in 0 or more than 1 exon)")
    if (last_exon_with_cds < first_exon_with_cds)
        stop("UCSC data anomaly in transcript ", tx_name,
             ": last exon with cds occurs before first exon with cds")
    exons_with_cds <- first_exon_with_cds:last_exon_with_cds
    cds_start0[exons_with_cds] <- exon_start0[exons_with_cds]
    cds_end1[exons_with_cds] <- exon_end1[exons_with_cds]
    cds_start0[first_exon_with_cds] <- cdsStart0
    cds_end1[last_exon_with_cds] <- cdsEnd1
    if (sum(cds_end1 - cds_start0, na.rm=TRUE) %% 3L != 0L)
        warning("UCSC data anomaly in transcript ", tx_name,
                ": the cds cumulative length is not a multiple of 3")
    list(cds_start0, cds_end1)
}

### 'exon_locs' must be the list of 2 lists returned by
### .extractExonLocsFromUCSCTxTable().
### Returns a named list with 2 elements. Each element is itself a list of
### integer vectors with eventually NAs. The 2 elements have the same
### "shape" as the elements in 'exon_locs' and the NAs occur at the same
### places in the 2 elements.
.extractCdsLocsFromUCSCTxTable <- function(ucsc_txtable, exon_locs)
{
    cdsStart <- ucsc_txtable$cdsStart
    cdsEnd <- ucsc_txtable$cdsEnd
    cds_start <- cds_end <- vector(mode="list", length=nrow(ucsc_txtable))
    for (i in seq_len(nrow(ucsc_txtable))) {
        startend <- .extractUCSCCdsStartEnd(cdsStart[i], cdsEnd[i],
                                            exon_locs$start[[i]],
                                            exon_locs$end[[i]],
                                            ucsc_txtable$name[i])
        cds_start[[i]] <- startend[[1L]]
        cds_end[[i]] <- startend[[2L]]
    }
    return(list(start=cds_start, end=cds_end))
}

.extractSplicingsFromUCSCTxTable <- function(ucsc_txtable, transcripts_tx_id)
{
    exon_count <- ucsc_txtable$exonCount
    splicings_tx_id <- rep.int(transcripts_tx_id, exon_count)
    if (length(exon_count) == 0L) {
        exon_rank <- exon_start <- exon_end <-
            cds_start <- cds_end <- integer(0)
    } else {
        if (min(exon_count) <= 0L)
            stop("UCSC data anomaly: 'ucsc_txtable$exonCount' contains ",
                 "non-positive values")
        exon_rank <- .makeExonRank(exon_count, ucsc_txtable$strand)
        exon_locs <- .extractExonLocsFromUCSCTxTable(ucsc_txtable)
        cds_locs <- .extractCdsLocsFromUCSCTxTable(ucsc_txtable, exon_locs)
        exon_start <- unlist(exon_locs$start) + 1L
        exon_end <- unlist(exon_locs$end)
        cds_start <- unlist(cds_locs$start) + 1L
        cds_end <- unlist(cds_locs$end)
    }
    splicings <- data.frame(
        tx_id=splicings_tx_id,
        exon_rank=exon_rank,
        exon_start=exon_start,
        exon_end=exon_end,
        cds_start=cds_start,
        cds_end=cds_end
    )
    splicings
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Preprocess the 'genes' data frame.
###

.makeUCSCGenes <- function(genes, ucsc_txtable)
{
    #genes <- genes[genes$tx_name %in% ucsc_txtable$name, ]
    genes
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'chrominfo' data frame.
###

.downloadChromInfoFromUCSC <- function(genome)
{
    url <- paste("http://hgdownload.cse.ucsc.edu/goldenPath/", genome,
                 "/database/chromInfo.txt.gz", sep="")
    destfile <- tempfile()
    download.file(url, destfile, quiet=TRUE)
    colnames <- c("chrom", "size", "fileName")
    ans <- read.table(destfile, sep="\t", quote="",
                      col.names=colnames, comment.char="")
    ans
}

.makeUCSCChrominfo <- function(genome)
{
    ucsc_chrominfotable <- .downloadChromInfoFromUCSC(genome)
    COL2CLASS <- c(
        chrom="character",
        size="integer"
    )
    ucsc_chrominfotable <- setDataFrameColClass(ucsc_chrominfotable, COL2CLASS,
                                                drop.extra.cols=TRUE)
    chrominfo <- data.frame(
        chrom=ucsc_chrominfotable$chrom,
        length=ucsc_chrominfotable$size
    )
    chrominfo
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Prepare the 'metadata' data frame.
###

.makeUCSCMetadata <- function(genome, tablename, gene_id_type, full_dataset)
{
    metadata <- data.frame(
        name=c("Data source", "Genome", "UCSC Table",
               "Type of Gene ID", "Full dataset"),
        value=c("UCSC", genome, tablename,
                gene_id_type, ifelse(full_dataset, "yes", "no"))
    )
    metadata
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTranscriptDbFromUCSC()
###

.makeTranscriptDbFromUCSCTxTable <- function(ucsc_txtable, genes,
                                             genome, tablename, gene_id_type,
                                             full_dataset)
{
    COL2CLASS <- c(
        name="character",
        chrom="factor",
        strand="factor",
        txStart="integer",
        txEnd="integer",
        cdsStart="integer",
        cdsEnd="integer",
        exonCount="integer",
        exonStarts="character",
        exonEnds="character"
    )
    ucsc_txtable <- setDataFrameColClass(ucsc_txtable, COL2CLASS,
                                         drop.extra.cols=TRUE)

    transcripts <- .extractTranscriptsFromUCSCTxTable(ucsc_txtable)
    splicings <- .extractSplicingsFromUCSCTxTable(ucsc_txtable,
                                                  transcripts$tx_id)
    genes <- .makeUCSCGenes(genes, ucsc_txtable)
    chrominfo <- .makeUCSCChrominfo(genome)
    metadata <- .makeUCSCMetadata(genome, tablename, gene_id_type, full_dataset)

    ## Call makeTranscriptDb().
    makeTranscriptDb(transcripts, splicings,
                     genes=genes, chrominfo=chrominfo, metadata=metadata)
}

### The 2 main tasks that makeTranscriptDbFromUCSC() performs are:
###   (1) download the data from UCSC into a data.frame (the getTable() call);
###   (2) store that data.frame in an SQLite db (the
###       .makeTranscriptDbFromUCSCTxTable() call).
### Speed:
###   - for genome="hg18" and tablename="knownGene":
###       (1) download takes about 40-50 sec.
###       (2) db creation takes about 30-35 sec.
makeTranscriptDbFromUCSC <- function(genome="hg18",
                                     tablename="knownGene",
                                     transcript_ids=NULL)
{
    if (!isSingleString(genome))
        stop("'genome' must be a single string")
    if (!isSingleString(tablename))
        stop("'tablename' must be a single string")
    track <- supportedUCSCtables()[tablename, "track"]
    if (is.na(track))
        stop("table \"", tablename, "\" is not supported")
    if (!is.null(transcript_ids)) {
        if (!is.character(transcript_ids) || any(is.na(transcript_ids)))
            stop("'transcript_ids' must be a character vector with no NAs")
    }
    session <- browserSession()
    genome(session) <- genome

    ## Download the transcript table.
    query <- ucscTableQuery(session, track, table=tablename,
                            names=transcript_ids)
    ucsc_txtable <- getTable(query)

    ## Download the tx_name-to-gene_id mapping.
    txname2geneid <- .fetchTxName2GeneIdMappingFromUCSC(session,
                                                        track, tablename)

    .makeTranscriptDbFromUCSCTxTable(ucsc_txtable, txname2geneid$genes,
                                     genome, tablename,
                                     txname2geneid$gene_id_type,
                                     full_dataset=is.null(transcript_ids))
}

