### =========================================================================
### makeTxDbFromUCSC()
### -------------------------------------------------------------------------


### makeTxDbFromUCSC() expects a UCSC transcript table to have at least
### the following columns:
.UCSC_TXCOL2CLASS <- c(
    name="character",
    chrom="character",
    strand="character",
    txStart="integer",
    txEnd="integer",
    cdsStart="integer",
    cdsEnd="integer",
    exonCount="integer",
    exonStarts="list",    # list of raw vectors
    exonEnds="list"       # list of raw vectors
)

### A thin wrapper around list_UCSC_primary_tables(genome, group="genes").
### Returns a data.frame with 3 columns ("tablename", "track", and
### "composite_track") and 1 row per primary table.
supportedUCSCtables <- function(genome="hg19",
                                url="https://genome.ucsc.edu/cgi-bin/")
{
    if (!missing(url))
        .Deprecated(msg="'url' argument is deprecated and was ignored")
    if (!isSingleString(genome))
        stop(wmsg("'genome' must be a single string"))
    df <- list_UCSC_primary_tables(genome, group="genes")
    expected_colnames <- c("primary_table", "track", "group", "composite_track")
    stopifnot(identical(colnames(df), expected_colnames))
    ans <- df[ , -3L]
    colnames(ans)[[1L]] <- "tablename"
    ans
}

### Can be used to quickly check that a combination of genome/tablename
### actually exists.
browseUCSCtrack <- function(genome="hg19",
                            tablename="knownGene",
                            url="https://genome.ucsc.edu/cgi-bin/")
{
    if (!isSingleString(genome))
        stop(wmsg("'genome' must be a single string"))
    if (!isSingleString(tablename))
        stop(wmsg("'tablename' must be a single string"))
    if (!isSingleString(url))
        stop(wmsg("'url' must be a single string"))
    url <- sprintf("%s/hgTrackUi?db=%s&g=%s", url, genome, tablename)
    ## Avoid "file association for 'https://...' not available or invalid"
    ## error during 'R CMD check' on some Windows systems (e.g. riesling1).
    if (interactive())
        browseURL(url)
}

.tablename2track <- function(tablename, genome)
{
    if (!isSingleString(tablename))
        stop(wmsg("'tablename' must be a single string"))
    supported_tables <- supportedUCSCtables(genome)
    idx <- match(tablename, supported_tables$tablename)
    if (is.na(idx))
        stop(wmsg("UCSC table \"", tablename, "\" is not supported"))
    supported_tables$track[idx]
}

### The table names above (unique key) must be used to name the top-level
### elements of the list below. If no suitable tx_name-to-gene_id mapping is
### available in the UCSC database for a supported table, then there is no
### entry in the list below for this table and makeTxDbFromUCSC() will leave
### the gene table empty.
.UCSC_TXNAME2GENEID_MAPDEFS <- list(
    knownGene=list(
        L2Rchain=list(
            c(tablename="knownToLocusLink",
              Lcolname="name",
              Rcolname="value")
        ),
        gene_id_type="Entrez Gene ID"
    ),
    refGene=list(
        L2Rchain=list(
            c(tablename="hgFixed.refLink",
              Lcolname="mrnaAcc",
              Rcolname="locusLinkId")
        ),
        gene_id_type="Entrez Gene ID"
    ),
    ncbiRefSeq=list(
        L2Rchain=list(
            c(tablename="ncbiRefSeqLink",
              Lcolname="id",
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
    ## UCSC changed its db schema in September 2011 and the ensGtp table became
    ## unavailable for some genomes (sacCer2, hg18, etc..., apparently for
    ## those assemblies that are not the latest). But the (new?) name2 column
    ## in the ensGene table seems to contain the Ensembl gene ids so the join
    ## with the ensGtp table is not needed anymore.
    #ensGene=list(
    #    L2Rchain=list(
    #        c(tablename="ensGtp",
    #          Lcolname="transcript",
    #          Rcolname="gene")
    #    ),
    #)
    ensGene=c(
        colname="name2",
        gene_id_type="Ensembl gene ID"
    ),
    lincRNAsTranscripts=c(
        colname="name",
        gene_id_type="Name of gene"
    ),
    wgEncodeGencodeManualV3=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeClassesV3",
              Lcolname="name",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeAutoV3=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeClassesV3",
              Lcolname="name",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePolyaV3=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeClassesV3",
              Lcolname="name",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeBasicV17=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV17",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeCompV17=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV17",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePseudoGeneV17=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV17",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencode2wayConsPseudoV17=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV17",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePolyaV17=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV17",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeBasicV14=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV14",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeCompV14=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV14",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePseudoGeneV14=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV14",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencode2wayConsPseudoV14=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV14",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePolyaV14=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV14",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeBasicV7=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV7",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodeCompV7=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV7",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePseudoGeneV7=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV7",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencode2wayConsPseudoV7=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV7",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    wgEncodeGencodePolyaV7=list(
        L2Rchain=list(
            c(tablename="wgEncodeGencodeAttrsV7",
              Lcolname="transcriptId",
              Rcolname="geneId")
        ),
        gene_id_type="Ensembl gene ID"
    ),
    flyBaseGene=list(
        L2Rchain=list(
            c(tablename="flyBaseIsoforms",
              Lcolname="transcript",
              Rcolname="clusterId"),
            c(tablename="flyBaseCanonical",
              Lcolname="clusterId",
              Rcolname="transcript")
        ),
        gene_id_type="Name of canonical transcript in cluster"
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
        gene_id_type="Name of canonical transcript in cluster"
    )
)

.get_txname2geneid_mapdef <- function(tablename)
    .UCSC_TXNAME2GENEID_MAPDEFS[[tablename]]

.fetch_UCSC_txtable <- function(genome, tablename, transcript_ids=NULL)
{
    if (is.null(transcript_ids)) {
        where <- NULL
    } else {
        where <- sprintf("name IN (%s)",
                         paste(paste0("'", transcript_ids, "'"), collapse=","))
    }
    columns <- names(.UCSC_TXCOL2CLASS)
    mapdef <- .get_txname2geneid_mapdef(tablename)
    if (is.character(mapdef))
        columns <- c(columns, mapdef[["colname"]])
    message("Download the ", tablename, " table ... ", appendLF=FALSE)
    ans <- UCSC_dbselect(genome, tablename, columns=columns, where=where)
    message("OK")
    ## DBI is returning blobs for exon starts and stops so the old check fails
    stopifnot(all(mapply(function(x, y) is(x, y), ans, .UCSC_TXCOL2CLASS)))
    ##current_classes <- head(sapply(ans, class),
    ##                        n=length(.UCSC_TXCOL2CLASS))
    ##stopifnot(identical(current_classes, .UCSC_TXCOL2CLASS))
    ans$exonStarts <- toListOfIntegerVectors(ans$exonStarts)
    ans$exonEnds <- toListOfIntegerVectors(ans$exonEnds)
    if (!identical(lengths(ans$exonStarts),
                   ans$exonCount))
        stop(wmsg("UCSC data anomaly in table ", genome, ".", tablename, ": ",
                  "columns exonStarts and exonCount are inconsistent"))
    if (!identical(lengths(ans$exonEnds),
                   ans$exonCount))
        stop(wmsg("UCSC data anomaly in table ", genome, ".", tablename, ": ",
                  "columns exonEnds and exonCount are inconsistent"))
    ans
}

.fetch_UCSC_table <- function(genome, tablename, columns=NULL)
{
    message("Download the ", tablename, " table ... ", appendLF=FALSE)
    if (tablename == "hgFixed.refLink") {
        genome <- "hgFixed"
        tablename <- "refLink"
    }
    ans <- UCSC_dbselect(genome, tablename, columns=columns)
    message("OK")
    ans
}

### The 2 functions below must return a named list with 2 elements:
###   $genes: data.frame with tx_name and gene_id cols;
###   $gene_id_type: single string.
.fetch_txname2geneid_from_UCSC <- function(genome, Ltablename, mapdef)
{
    nlink <- length(mapdef$L2Rchain)
    for (i in seq_len(nlink)) {
        L2Rlink <- mapdef$L2Rchain[[i]]
        tablename <- L2Rlink[["tablename"]]
        Lcolname <- L2Rlink[["Lcolname"]]
        Rcolname <- L2Rlink[["Rcolname"]]
        ## The tables involved in the "left join" don't necessarily belong
        ## to the track of the leftmost table (e.g.
        ## "wgEncodeGencodeAttrsV17" table does NOT belong to the "GENCODE
        ## Genes V17" track) so we don't need the track information to fetch
        ## these tables.
        ucsc_table <- .fetch_UCSC_table(genome, tablename,
                                        columns=c(Lcolname, Rcolname))
        if (!all(has_col(ucsc_table, c(Lcolname, Rcolname))))
            stop(wmsg("expected cols \"", Lcolname, "\" or/and \"",
                      Rcolname, "\" not found in table ", tablename))
        Lcol <- ucsc_table[[Lcolname]]
        Rcol <- ucsc_table[[Rcolname]]
        if (!is.character(Lcol))
            Lcol <- as.character(Lcol)
        if (!is.character(Rcol))
            Rcol <- as.character(Rcol)
        if (i == 1L) {
            tmp <- data.frame(Lcol=Lcol, Rcol=Rcol, stringsAsFactors=FALSE)
        } else {
            name2val <- Rcol
            names(name2val) <- Lcol
            tmp <- joinDataFrameWithName2Val(tmp, "Rcol", name2val, "Rcol")
        }
    }
    genes <- data.frame(tx_name=tmp$Lcol,
                        gene_id=tmp$Rcol,
                        stringsAsFactors=FALSE)
    gene_id_type <- mapdef$gene_id_type
    list(genes=genes, gene_id_type=gene_id_type)
}

.extract_txname2geneid_from_UCSC_txtable <- function(ucsc_txtable, mapdef)
{
    genes <- data.frame(tx_name=ucsc_txtable[["name"]],
                        gene_id=ucsc_txtable[[mapdef[["colname"]]]],
                        stringsAsFactors=FALSE)
    gene_id_type <- mapdef[["gene_id_type"]]
    list(genes=genes, gene_id_type=gene_id_type)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the 'transcripts' data frame from UCSC table.
###

.extract_transcripts_from_UCSC_txtable <- function(ucsc_txtable)
{
    message("Extract the 'transcripts' data frame ... ", appendLF=FALSE)
    tx_id <- seq_len(nrow(ucsc_txtable))
    tx_name <- ucsc_txtable$name
    tx_chrom <- ucsc_txtable$chrom
    tx_strand <- ucsc_txtable$strand
    tx_start <- ucsc_txtable$txStart + 1L
    tx_end <- ucsc_txtable$txEnd
    ans <- data.frame(
        tx_id=tx_id,
        tx_name=tx_name,
        tx_chrom=tx_chrom,
        tx_strand=tx_strand,
        tx_start=tx_start,
        tx_end=tx_end,
        stringsAsFactors=FALSE
    )
    message("OK")
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the 'splicings' data frame from UCSC table.
###

### 'cdsStart0', 'cdsEnd1': single integers (resp. 0-based and 1-based).
### 'exon_start0', 'exon_end1': integer vectors of equal lengths (resp.
### 0-based and 1-based) and with no NAs.
### Returns a list with 2 elements, each of them being an integer vector of
### the same length as 'exon_start0' (or 'exon_end1') that may contain NAs.
### Notes (using genome="hg18"):
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
###   (3) All transcripts in the knownGene and ccdsGene tables have a cds
###       cumulative length that is a multiple of 3 (the former doesn't even
###       have the cdsStartStat/cdsStartEnd columns). For hg19 ccdsGene, this
###       is not true anymore :-/
### TODO: Investigate (1) and (2).
.extract_UCSC_cds_start_end <- function(cdsStart0, cdsEnd1,
                                        exon_start0, exon_end1, tx_name)
{
    cds_start0 <- cds_end1 <- integer(length(exon_start0))
    cds_start0[] <- NA_integer_
    cds_end1[] <- NA_integer_
    if (cdsStart0 >= cdsEnd1)
        return(list(cds_start0, cds_end1, FALSE))
    first_exon_with_cds <- which(exon_start0 <= cdsStart0
                                 & cdsStart0 < exon_end1)
    if (length(first_exon_with_cds) != 1L)
        stop(wmsg("UCSC data ambiguity in transcript ", tx_name,
                  ": cannot determine first exon with cds ('cdsStart' ",
                  "falls in 0 or more than 1 exon)"))
    last_exon_with_cds <- which(exon_start0 < cdsEnd1
                                & cdsEnd1 <= exon_end1)
    if (length(last_exon_with_cds) != 1L)
        stop(wmsg("UCSC data ambiguity in transcript ", tx_name,
                  ": cannot determine last exon with cds ('cdsEnd' ",
                  "falls in 0 or more than 1 exon)"))
    if (last_exon_with_cds < first_exon_with_cds)
        stop(wmsg("UCSC data anomaly in transcript ", tx_name,
                  ": last exon with cds occurs before first exon with cds"))
    exons_with_cds <- first_exon_with_cds:last_exon_with_cds
    cds_start0[exons_with_cds] <- exon_start0[exons_with_cds]
    cds_end1[exons_with_cds] <- exon_end1[exons_with_cds]
    cds_start0[first_exon_with_cds] <- cdsStart0
    cds_end1[last_exon_with_cds] <- cdsEnd1
    ## NAs are OK in here since they indicate the absence of any CDS
    ## (which is common and nothing to write home about)
    ## changed from 50K to 19.3K ...   but the 'bad' ones are not present?
    bad <- sum(cds_end1 - cds_start0, na.rm=TRUE) %% 3L != 0L
    list(cds_start0, cds_end1, bad)
}

### Return a named list with 2 list elements, each of which is itself a list
### of integer vectors with eventually NAs. The 2 elements have the same
### "shape" as ucsc_txtable$exonStarts and ucsc_txtable$exonEnds and the NAs
### in them occur at the same places in the 2 elements.
.extract_cds_locs_from_UCSC_txtable <- function(ucsc_txtable)
{
    start_end <- Map(.extract_UCSC_cds_start_end,
                     ucsc_txtable$cdsStart, ucsc_txtable$cdsEnd,
                     ucsc_txtable$exonStarts, ucsc_txtable$exonEnds,
                     ucsc_txtable$name)
    bad <- sapply(start_end, "[[", 3L)
    if (any(bad)) {
        bad_cds <- ucsc_txtable$name[bad]
        msg <- sprintf("UCSC data anomaly in %d transcript(s):
            the cds cumulative length is not a multiple of 3
            for transcripts %s", length(bad_cds),
            paste(sQuote(bad_cds), collapse=" "))
        warning(wmsg(paste(strwrap(msg, exdent=2L), collapse="\n")))
    }
    list(start=sapply(start_end, "[[", 1L),
         end=sapply(start_end, "[[", 2L))
}

.extract_splicings_from_UCSC_txtable <- function(ucsc_txtable,
                                                 transcripts_tx_id)
{
    message("Extract the 'splicings' data frame ... ", appendLF=FALSE)
    exon_count <- ucsc_txtable$exonCount
    splicings_tx_id <- rep.int(transcripts_tx_id, exon_count)
    if (length(exon_count) == 0L) {
        exon_rank <- exon_start <- exon_end <-
            cds_start <- cds_end <- integer(0)
    } else {
        if (min(exon_count) <= 0L)
            stop(wmsg("UCSC data anomaly: 'ucsc_txtable$exonCount' contains ",
                      "non-positive values"))
        exon_rank <- make_exon_rank_col(exon_count, ucsc_txtable$strand)
        cds_locs <- .extract_cds_locs_from_UCSC_txtable(ucsc_txtable)
        exon_start <- unlist(ucsc_txtable$exonStarts) + 1L
        exon_end <- unlist(ucsc_txtable$exonEnds)
        cds_start <- unlist(cds_locs$start) + 1L
        cds_end <- unlist(cds_locs$end)
    }
    ans <- data.frame(
        tx_id=splicings_tx_id,
        exon_rank=exon_rank,
        exon_start=exon_start,
        exon_end=exon_end,
        cds_start=cds_start,
        cds_end=cds_end,
        stringsAsFactors=FALSE
    )
    message("OK")
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Preprocess the 'genes' data frame.
###

.make_UCSC_genes <- function(genes, ucsc_txtable)
{
    #genes <- S4Vectors:::extract_data_frame_rows(genes,
    #                             genes$tx_name %in% ucsc_txtable$name)
    genes
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Download and preprocess the 'chrominfo' data frame.
###

.make_UCSC_chrominfo <- function(genome, circ_seqs=NULL,
        goldenPath.url=getOption("UCSC.goldenPath.url"))
{
    message("Download and preprocess the 'chrominfo' data frame ... ",
            appendLF=FALSE)
    warning_tip1 <- paste0("You can see this by calling isCircular() on the ",
                           "TxDb object returned by makeTxDbFromUCSC().")
    warning_tip2 <- paste0("See '?makeTxDbFromUCSC' in the GenomicFeatures ",
                           "package for more information.")
    chrominfo <- get_and_fix_chrom_info_from_UCSC(genome,
                                   goldenPath.url=goldenPath.url,
                                   circ_seqs=circ_seqs,
                                   warning_tip1=warning_tip1,
                                   warning_tip2=warning_tip2)
    ## Prepare the data frame that will be passed to the 'chrominfo' argument
    ## of makeTxDb(). Note that the naming convention for the columns in this
    ## data frame differs from what's used in the data frame returned by
    ## get_and_fix_chrom_info_from_UCSC().
    ans <- data.frame(
        chrom=chrominfo[ , "chrom"],
        length=chrominfo[ , "size"],
        is_circular=chrominfo[ , "circular"],
        stringsAsFactors=FALSE
    )
    message("OK")
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Prepare the 'metadata' data frame.
###

.prepare_UCSC_metadata <- function(genome, tablename, track, gene_id_type,
                                   full_dataset,
                                   taxonomyId=NA, miRBaseBuild=NA)
{
    message("Prepare the 'metadata' data frame ... ", appendLF=FALSE)
    if (!isSingleStringOrNA(miRBaseBuild))
        stop(wmsg("'miRBaseBuild' must be a a single string or NA"))
    organism <- lookup_organism_by_UCSC_genome(genome)
    if (is.na(taxonomyId)) {
        taxonomyId <- GenomeInfoDb:::lookup_tax_id_by_organism(organism)
    } else {
        GenomeInfoDb:::check_tax_id(taxonomyId)
    }

    ans <- data.frame(
        name=c("Data source", "Genome", "Organism", "Taxonomy ID",
               "UCSC Table", "UCSC Track",
               "Resource URL", "Type of Gene ID",
               "Full dataset",
               "miRBase build ID"),
        value=c("UCSC", genome, organism, taxonomyId,
                tablename, track,
                "https://genome.ucsc.edu/", gene_id_type,
                ifelse(full_dataset, "yes", "no"),
                miRBaseBuild),
        stringsAsFactors=FALSE
    )
    message("OK")
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTxDbFromUCSC()
###

.make_TxDb_from_UCSC_txtable <- function(ucsc_txtable, genes,
        genome, tablename, track, gene_id_type,
        full_dataset,
        circ_seqs=NULL,
        goldenPath.url=getOption("UCSC.goldenPath.url"),
        taxonomyId=NA,
        miRBaseBuild=NA)
{
    strand_is_dot <- ucsc_txtable$strand == "."
    if (any(strand_is_dot)) {
        msg <- sprintf("dropped %d transcript(s) for which strand
            was not set (i.e. was set to '.')", sum(strand_is_dot))
        warning(wmsg(msg))
        keep_idx <- which(!strand_is_dot)
        ucsc_txtable <- S4Vectors:::extract_data_frame_rows(ucsc_txtable,
                                                            keep_idx)
    }

    transcripts <- .extract_transcripts_from_UCSC_txtable(ucsc_txtable)
    splicings <- .extract_splicings_from_UCSC_txtable(ucsc_txtable,
                                                      transcripts$tx_id)
    genes <- .make_UCSC_genes(genes, ucsc_txtable)
    chrominfo <- .make_UCSC_chrominfo(genome, circ_seqs, goldenPath.url)
    metadata <- .prepare_UCSC_metadata(genome, tablename, track, gene_id_type,
                                       full_dataset,
                                       taxonomyId,  miRBaseBuild)
    ## Jan 2019 -- The refGene tables in the hg19 and hg38 UCSC databases were
    ## last updated in Nov 2018 and now contain transcripts located on
    ## sequences that don't belong to the underlying genomes (GRCh37 and GRCh38
    ## respectively). More precisely some transcripts in these tables now
    ## belong to patched versions of these genomes: GRCh37.p13 for hg19 and
    ## GRCh38.p11 for hg38. This causes the makeTxDbFromUCSC() errors reported
    ## here:
    ##   https://github.com/Bioconductor/GenomicFeatures/issues/14
    ##   https://support.bioconductor.org/p/117265/
    ##   https://support.bioconductor.org/p/114901/
    ## The current fix is to drop these foreign transcripts with a warning.
    if (genome %in% c("hg19", "hg38") && tablename == "refGene")
        on.foreign.transcripts <- "drop"
    else
        on.foreign.transcripts <- "error"

    message("Make the TxDb object ... ", appendLF=FALSE)
    ans <- makeTxDb(transcripts, splicings, genes=genes,
                    chrominfo=chrominfo, metadata=metadata,
                    reassign.ids=TRUE,
                    on.foreign.transcripts=on.foreign.transcripts)
    message("OK")
    ans
}

### Some timings (as of Jan 31, 2018, GenomicFeatures 1.31.7):
###          |             |    nb of    |
###   genome |   tablename | transcripts | time (s)
###   ---------------------------------------------
###     hg18 |   knownGene |       66803 |     37.2
###     hg18 |     refGene |       68178 |     42.6
###     hg19 |   knownGene |       82960 |     44.9
###     hg19 |     refGene |       69998 |     45.5
###     hg38 |     ensGene |      204940 |     63.7
###     hg19 |    ccdsGene |       28856 |     29.2
###     hg19 | xenoRefGene |      177746 |    114.1
###     hg38 |   knownGene |      197782 |     53.9
###     hg38 |     refGene |       74673 |     38.4
###      dm3 | flyBaseGene |       21236 |     28.9
###  sacCer2 |     sgdGene |        6717 |     22.6
###  sacCer3 |     ensGene |        7126 |     18.1
makeTxDbFromUCSC <- function(genome="hg19",
        tablename="knownGene",
        transcript_ids=NULL,
        circ_seqs=NULL,
        url="https://genome.ucsc.edu/cgi-bin/",
        goldenPath.url=getOption("UCSC.goldenPath.url"),
        taxonomyId=NA,
        miRBaseBuild=NA)
{
    if (!requireNamespace("RMariaDB", quietly=TRUE))
        stop(wmsg("Couldn't load the RMariaDB package. ",
                  "You need to install the RMariaDB package ",
                  "in order to use makeTxDbFromUCSC()."))

    if (!is.null(transcript_ids)) {
        if (!is.character(transcript_ids) || any(is.na(transcript_ids)))
            stop(wmsg("'transcript_ids' must be a ",
                      "character vector with no NAs"))
    }
    if (!missing(url))
        .Deprecated(msg="'url' argument is deprecated and was ignored")
    if (!isSingleString(goldenPath.url))
        stop(wmsg("'goldenPath.url' must be a single string"))

    track <- .tablename2track(tablename, genome)

    ## Download the transcript table.
    ucsc_txtable <- .fetch_UCSC_txtable(genome, tablename,
                                        transcript_ids=transcript_ids)

    ## Get the tx_name-to-gene_id mapping.
    mapdef <- .get_txname2geneid_mapdef(tablename)
    if (is.null(mapdef)) {
        txname2geneid <- list(genes=NULL, gene_id_type="no gene ids")
    } else if (is.list(mapdef)) {
        txname2geneid <- .fetch_txname2geneid_from_UCSC(
                                 genome,
                                 tablename, mapdef)
    } else if (is.character(mapdef)) {
        txname2geneid <- .extract_txname2geneid_from_UCSC_txtable(
                                 ucsc_txtable, mapdef)
    } else {
        stop(wmsg("GenomicFeatures internal error: invalid 'mapdef'"))
    }
    .make_TxDb_from_UCSC_txtable(ucsc_txtable, txname2geneid$genes,
                                 genome, tablename, track,
                                 txname2geneid$gene_id_type,
                                 full_dataset=is.null(transcript_ids),
                                 circ_seqs=circ_seqs,
                                 goldenPath.url=goldenPath.url,
                                 taxonomyId=taxonomyId,
                                 miRBaseBuild=miRBaseBuild)
}

