### =========================================================================
### Making TranscriptDb objects
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Saving/loading.
###

saveFeatures <- function(x, file)
{
    if (!is(x, "TranscriptDb"))
        stop("'x' must be a TranscriptDb object")
    if (!isSingleString(file))
        stop("'file' must be a single string")
    sqliteCopyDatabase(x@conn, file)
}

### FIXME: loadFeatures() needs to put the db back into memory (it's currently
### returning a TranscriptDb object that points to the on-disk db).
loadFeatures <- function(file)
{
    if (!isSingleString(file))
        stop("'file' must be a single string")
    if(!file.exists(file))
        stop("file '", file, "' does not exist")
    conn <- dbConnect(SQLite(), file)
    new("TranscriptDb", conn=conn)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions for .makeTranscriptDb() / makeTranscriptDb().
###

.argAsCharacterFactorWithNoNAs <- function(arg, argname)
{
    if (is.factor(arg)) {
        if (is.character(levels(arg)) && !any(is.na(arg)))
            return(arg)
        stop("'", argname, "' must be a character vector/factor with no NAs")
    }
    if (!is.character(arg) || any(is.na(arg)))
        stop("'", argname, "' must be a character vector/factor with no NAs")
    as.factor(arg)
}

.argAsIntegerWithNoNAs <- function(arg, argname)
{
    if (!is.numeric(arg) || any(is.na(arg)))
        stop("'", argname, "' must be an integer vector with no NAs")
    if (!is.integer(arg))
        arg <- as.integer(arg)
    arg
}

.makeInternalIdsFromExternalIds <- function(external_id)
    as.integer(factor(external_id))

.makeInternalIdsForUniqueLocs <- function(chrom, strand, start, end)
{
    x <- data.frame(chrom, strand, start, end, stringsAsFactors=FALSE)
    makeIdsForUniqueDataFrameRows(x)
}

### Because we use SQLite "rtree" feature, .writeFeatureCoreTables() creates
### the 5 following tables:
###   (1) '<feature>s'
###   (2) '<feature>s_rtree'
###   (3) '<feature>s_rtree_node'
###   (4) '<feature>s_rtree_parent'
###   (5) '<feature>s_rtree_rowid'
### Note that only (1) and (2) are explicitely created. (3), (4) and (5) are
### automatically created by the SQLite engine.
.writeFeatureCoreTables <- function(conn,
                                    feature,
                                    internal_id,
                                    external_id,
                                    name,
                                    chrom,
                                    strand,
                                    start,
                                    end,
                                    colnames)
{
    if (is.null(name))
        name <- rep.int(NA_character_, length(internal_id))
    table <- data.frame(
        internal_id=internal_id,
        external_id=external_id,
        name=name,
        chrom=chrom,
        strand=strand,
        start=start,
        end=end,
        stringsAsFactors=FALSE)
    table <- unique(table)
    ## Create the '<feature>s' table.
    sql <- c(
        "CREATE TABLE ", feature, "s (\n",
        "  ", colnames[1L], " INTEGER PRIMARY KEY,\n",
        "  ", colnames[2L], " TEXT UNIQUE NOT NULL,\n",
        "  ", colnames[3L], " TEXT NULL,\n",
        "  ", colnames[4L], " TEXT NOT NULL,\n",
        "  ", colnames[5L], " TEXT NOT NULL\n",
        ")")
    res <- dbSendQuery(conn, paste(sql, collapse=""))
    dbClearResult(res)
    ## Fill the '<feature>s' table.
    ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
    ## when the nb of rows to insert is 0, hence the following test.
    if (nrow(table) != 0L) {
        sql <- c("INSERT INTO ", feature, "s VALUES (?,?,?,?,?)")
        dbBeginTransaction(conn)
        res <- dbSendPreparedQuery(conn, paste(sql, collapse=""),
                                   table[1:5])
        dbClearResult(res)
        dbCommit(conn)
    }
    ## Create the '<feature>s_rtree' table.
    sql <- c(
        "CREATE VIRTUAL TABLE ", feature, "s_rtree USING rtree (\n",
        "  ", colnames[1L], " INTEGER PRIMARY KEY,\n",
        "  ", colnames[6L], " INTEGER, -- NOT NULL is implicit in rtree\n",
        "  ", colnames[7L], " INTEGER  -- NOT NULL is implicit in rtree\n",
        "  -- FOREIGN KEY (", colnames[1L], ") REFERENCES ", feature, "s\n",
        ")")
    res <- dbSendQuery(conn, paste(sql, collapse=""))
    dbClearResult(res)
    ## Fill the '<feature>s_rtree' table.
    if (nrow(table) != 0L) {
        sql <- c("INSERT INTO ", feature, "s_rtree VALUES (?,?,?)")
        dbBeginTransaction(conn)
        res <- dbSendPreparedQuery(conn, paste(sql, collapse=""),
                                   table[c(1L, 6:7)])
        dbClearResult(res)
        dbCommit(conn)
    }
}

.writeSplicingsTable <- function(conn,
                                 internal_tx_id,
                                 exonRank,
                                 internal_exon_id,
                                 internal_CDS_id)
{
    table <- data.frame(
        internal_tx_id=internal_tx_id,
        exonRank=exonRank,
        internal_exon_id=internal_exon_id,
        internal_CDS_id=internal_CDS_id,
        stringsAsFactors=FALSE)
    table <- unique(table)
    ## Create the 'splicings' table.
    sql <- c(
        "CREATE TABLE splicings (\n",
        "  _tx_id INTEGER NOT NULL,\n",
        "  exon_rank INTEGER NOT NULL,\n",
        "  _exon_id INTEGER NOT NULL,\n",
        "  _CDS_id INTEGER NULL,\n",
        "  UNIQUE (_tx_id, exon_rank),\n",
        "  FOREIGN KEY (_tx_id) REFERENCES transcripts,\n",
        "  FOREIGN KEY (_exon_id) REFERENCES exons,\n",
        "  FOREIGN KEY (_CDS_id) REFERENCES CDSs\n",
        ")")
    res <- dbSendQuery(conn, paste(sql, collapse=""))
    dbClearResult(res)
    ## Fill the 'splicings' table.
    ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
    ## when the nb of rows to insert is 0, hence the following test.
    if (nrow(table) != 0L) {
        sql <- "INSERT INTO splicings VALUES (?,?,?,?)"
        dbBeginTransaction(conn)
        res <- dbSendPreparedQuery(conn, sql, table)
        dbClearResult(res)
        dbCommit(conn)
    }
}

.writeGenesTable <- function(conn, geneId, internal_tx_id)
{
    table <- data.frame(
        geneId=geneId,
        internal_tx_id=internal_tx_id,
        stringsAsFactors=FALSE)
    table <- unique(table)
    table <- table[!is.na(table$geneId), ]
    ## Create the 'genes' table.
    sql <- c(
        "CREATE TABLE genes (\n",
        "  gene_id TEXT NOT NULL,\n",
        "  _tx_id INTEGER NOT NULL,\n",
        "  UNIQUE (gene_id, _tx_id),\n",
        "  FOREIGN KEY (_tx_id) REFERENCES transcripts\n",
        ")")
    res <- dbSendQuery(conn, paste(sql, collapse=""))
    dbClearResult(res)
    ## Fill the 'genes' table.
    ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
    ## when the nb of rows to insert is 0, hence the following test.
    if (nrow(table) != 0L) {
        sql <- "INSERT INTO genes VALUES (?,?)"
        dbBeginTransaction(conn)
        res <- dbSendPreparedQuery(conn, sql, table)
        dbClearResult(res)
        dbCommit(conn)
    }
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### .makeTranscriptDb().
###

### WORK-IN-PROGRESS
.makeTranscriptDb <- function(transcript, exon, cds, ...)
{
    stop("WORK-IN-PROGRESS")
    args <- list(transcript, exon, cds, ...)
    if (!all(sapply(args, is.data.frame)))
        stop("all args must be data frames")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTranscriptDb().
###

### Creates a TranscriptDb instance from vectors of data.
### All vectors must be of equal length.  The i-th element of the
### vectors represent data for a single exon. Transcript data will
### be repeated over all exons within the trancsript.
makeTranscriptDb <- function(geneId,
                        txId, txName=NULL, txChrom, txStrand, txStart, txEnd,
                        cdsStart, cdsEnd,
                        exonStart, exonEnd, exonRank, exonId=NULL)
{
    if (!is.character(geneId))
        stop("'geneId' must be a character vector")
    if (!is.character(txId) || any(is.na(txId)))
        stop("'txId' must be a character vector with no NAs")
    if (!is.null(txName) && !is.character(txName))
        stop("when supplied, 'txName' must be a character vector")
    txChrom <- .argAsCharacterFactorWithNoNAs(txChrom, "txChrom")
    txStrand <- .argAsCharacterFactorWithNoNAs(txStrand, "txStrand")
    txStart <- .argAsIntegerWithNoNAs(txStart, "txStart")
    txEnd <- .argAsIntegerWithNoNAs(txEnd, "txEnd")
    #cdsStart <- .argAsIntegerWithNoNAs(cdsStart, "cdsStart")
    #cdsEnd <- .argAsIntegerWithNoNAs(cdsEnd, "cdsEnd")
    exonStart <- .argAsIntegerWithNoNAs(exonStart, "exonStart")
    exonEnd <- .argAsIntegerWithNoNAs(exonEnd, "exonEnd")
    exonRank <- .argAsIntegerWithNoNAs(exonRank, "exonRank")
    internal_tx_id <- .makeInternalIdsFromExternalIds(txId) 
    if (is.null(exonId)) {
        internal_exon_id <- .makeInternalIdsForUniqueLocs(
                                txChrom, txStrand, exonStart, exonEnd)
        exonId <- as.character(internal_exon_id)
    } else {
        if (!is.character(exonId) || any(is.na(exonId)))
            stop("when supplied, 'exonId' must be a character vector ",
                 "with no NAs")
        internal_exon_id <- .makeInternalIdsFromExternalIds(exonId)
    }

    conn <- dbConnect(SQLite(), dbname="") ## we'll write the db to a temp file
    .writeFeatureCoreTables(conn, "transcript",
        internal_tx_id, txId, txName, txChrom, txStrand, txStart, txEnd,
        c("_tx_id", "tx_id", "tx_name", "tx_chrom", "tx_strand",
          "tx_start", "tx_end"))
    .writeFeatureCoreTables(conn, "exon",
        internal_exon_id, exonId, NULL, txChrom, txStrand, exonStart, exonEnd,
        c("_exon_id", "exon_id", "exon_name", "exon_chrom", "exon_strand",
          "exon_start", "exon_end"))
    ## CDS table is temporarily left empty.
    .writeFeatureCoreTables(conn, "CDS",
        integer(0), character(0), NULL,
        character(0), character(0), integer(0), integer(0),
        c("_CDS_id", "CDS_id", "CDS_name", "CDS_chrom", "CDS_strand",
          "CDS_start", "CDS_end"))
    internal_CDS_id <- rep.int(NA_integer_, length(internal_tx_id))
    .writeSplicingsTable(conn,
        internal_tx_id, exonRank, internal_exon_id, internal_CDS_id)
    .writeGenesTable(conn, geneId, internal_tx_id)
    new("TranscriptDb", conn=conn)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTranscriptDbFromUCSC() and UCSC-specific helper functions.
###

### Could belong to rtracklayer.
.getOrganismFromUCSCgenome <- function(genome)
{
    ## Work around a bug in ucscGenomes().
    if (genome == "hg19")
        return("Human")
    ucsc_genomes <- ucscGenomes()
    which_genome <- match(genome, as.character(ucsc_genomes[,"db"]))
    if (is.na(which_genome))
        stop("'genome' not in ucscGenomes()")
    as.character(ucsc_genomes[which_genome, "organism"])
}

### TODO: This mapping from organism to org package could be generically
### useful in annotate.
.getOrgPkgFromOrganism <- function(organism)
{
    ### TODO: This mapping from organism to orgcode could be generically
    ### useful in annotate.
    ORGANISM2ORGCODE <- c(
        `A. gambiae`="Ag",
        `Cow`="Bt",
        `C. elegans`="Ce",
        `Dog`="Cf",
        `D. melanogaster`="Dm",
        `Zebrafish`="Dr",
        ## "EcK12"
        ## "EcSakai"
        `Chicken`="Gg",
        `Human`="Hs",
        `Mouse`="Mm",
        `Rhesus`="Mmu",
        `Chimp`="Pt",
        `Rat`="Rn"
        ## "Ss"
        ## "Xl"
    )
    orgcode <- unname(ORGANISM2ORGCODE[organism]) 
    if (is.na(orgcode))
        stop("no org.*.eg.db package for ", organism)
    paste("org", orgcode, "eg.db", sep=".")
}

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

.makeTranscriptDbFromUCSCTxTable <- function(ucsc_txtable, organism, track)
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
    ## Map the transcript IDs in 'ucsc_txtable$name' to the corresponding
    ## Entrez Gene IDs. Depending on the value of 'track' ("knownGene",
    ## "refGene" or "ensGene") this is done by using either the UCSCKG, the
    ## REFSEQ or the ENSEMBLTRANS map from the appropriate org package.
    orgpkg <- .getOrgPkgFromOrganism(organism)
    TRACK2MAPNAME <- c(
        knownGene="UCSCKG",
        refGene="REFSEQ",
        ensGene="ENSEMBLTRANS"
    )
    mapname <- unname(TRACK2MAPNAME[track])
    if (is.na(mapname))
        stop("don't know which map in ", orgpkg, " to use to map the ",
             "transcript IDs in track \"", track, "\" to Entrez Gene IDs")
    map <- suppressMessages(getAnnMap(mapname, orgpkg))
    ## Here we make the assumption that each transcript ID is mapped to at
    ## most one Entrez Gene ID, which will probably be true most of the
    ## time. For now we raise an error if this is not the case, but we might
    ## also want to handle this nicely in the future.
    tx2EG <- mget(ucsc_txtable$name, revmap(map), ifnotfound=NA_character_)
    if (any(sapply(tx2EG, length) > 1L))
        stop("some transcript IDs are mapped to multiple Entrez Gene IDs")
    EGs <- unname(unlist(tx2EG))  # can contain NAs
    ## Exon starts and ends are multi-valued fields (comma-separated) that
    ## need to be expanded.
    exon_count <- ucsc_txtable$exonCount
    if (min(exon_count) <= 0L)
        stop("'ucsc_txtable$exonCount' contains non-positive values")
    exonStarts <- strsplitAsListOfIntegerVectors(ucsc_txtable$exonStarts)
    if (!identical(elementLengths(exonStarts), exon_count))
        stop("'ucsc_txtable$exonStarts' inconsistent ",
             "with 'ucsc_txtable$exonCount'")
    exonEnds <- strsplitAsListOfIntegerVectors(ucsc_txtable$exonEnds)
    if (!identical(elementLengths(exonEnds), exon_count))
        stop("'ucsc_txtable$exonEnds' inconsistent ",
             "with 'ucsc_txtable$exonCount'")
    txName <- rep.int(ucsc_txtable$name, exon_count)
    ## For some tracks (e.g. knowGene), the 'name' col in the UCSC db is
    ## known not be a unique transcript identifier. But that's not always
    ## the case! For example, the refGene track uses the same transcript
    ## name for different transcripts. In that case, we need to generate
    ## our own transcript ids.
    ## TODO: Is it safe to treat the ensGene track like the knownGene track?
    if (track == "knownGene") {
        txId <- txName
    } else {
        txId <- as.character(rep.int(seq_len(nrow(ucsc_txtable)), exon_count))
    }
    txtable <- data.frame(
        geneId=rep.int(EGs, exon_count),
        txId=txId,
        txName=txName,
        txChrom=rep.int(ucsc_txtable$chrom, exon_count),
        txStrand=rep.int(ucsc_txtable$strand, exon_count),
        txStart=rep.int(ucsc_txtable$txStart, exon_count) + 1L,
        txEnd=rep.int(ucsc_txtable$txEnd, exon_count),
        cdsStart=rep.int(ucsc_txtable$cdsStart, exon_count) + 1L,
        cdsEnd=rep.int(ucsc_txtable$cdsEnd, exon_count),
        exonStart=unlist(exonStarts) + 1L,
        exonEnd=unlist(exonEnds),
        exonRank=.makeExonRank(exon_count, ucsc_txtable$strand),
        stringsAsFactors=FALSE)
    do.call(makeTranscriptDb, txtable)
}

### The 2 main tasks that makeTranscriptDbFromUCSC() performs are:
###   (1) download the data from UCSC into a data.frame (the getTable() call);
###   (2) store that data.frame in an SQLite db (the
###       .makeTranscriptDbFromUCSCTxTable() call).
### Speed:
###   - for genome="hg18" and track="knownGene":
###       (1) download takes about 40-50 sec.
###       (2) db creation takes about 30-35 sec.
### TODO: Support for track="refGene" is currently disabled because this
### track doesn't seem to contain anything that could be used as a unique
### transcript ID. So we need to either find a way to retrieve (or generate)
### this unique ID from somewhere else or drop "refGene" from the list of
### valid values for the 'track' arg.
### FIXME: Support for track="ensGene" is currently broken because some
### transcript IDs are mapped to multiple Entrez Gene IDs.
makeTranscriptDbFromUCSC <- function(genome="hg18",
                                     track=c("knownGene","refGene","ensGene"))
{
    if (!isSingleString(genome))
        stop("'genome' must be a single string")
    organism <- .getOrganismFromUCSCgenome(genome)
    track <- match.arg(track)
    if (track == "refGene")
        stop("track \"refGene\" is not currently supported, sorry!")
    if (track == "knownGene" && !(organism %in% c("Human", "Mouse", "Rat")))
        stop("UCSC \"knownGene\" track is only supported ",
             "for Human, Mouse and Rat")
    session <- browserSession()
    genome(session) <- genome
    query <- ucscTableQuery(session, track)
    ucsc_txtable <- getTable(query)  # download the data
    .makeTranscriptDbFromUCSCTxTable(ucsc_txtable, organism, track)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTranscriptDbFromBiomart().
###
### For people who want to tap biomaRt.
### Typical use:
###   txdb <- makeTranscriptDbFromBiomart(biomart="ensembl",
###                                       dataset="hsapiens_gene_ensembl")
### Speed:
###   - for biomart="ensembl" and dataset="hsapiens_gene_ensembl":
###       (1) download takes about 8 min.
###       (2) db creation takes about 60-65 sec.
###

makeTranscriptDbFromBiomart <- function(biomart="ensembl",
                                        dataset="hsapiens_gene_ensembl",
                                        ensembl_transcript_ids=NULL)
{
    mart <- useMart(biomart=biomart, dataset=dataset)
    attributes <- c("ensembl_gene_id",
                    "ensembl_transcript_id",
                    "chromosome_name",
                    "strand",
                    "transcript_start",
                    "transcript_end",
                    "cds_start",
                    "cds_end",
                    "ensembl_exon_id",
                    "exon_chrom_start",
                    "exon_chrom_end",
                    "rank")
    if (is.null(ensembl_transcript_ids)) {
        bm_txtable <- getBM(attributes, mart=mart)
    } else {
        if (!is.character(ensembl_transcript_ids)
         || any(is.na(ensembl_transcript_ids)))
            stop("'ensembl_transcript_ids' must be ",
                 "a character vector with no NAs")
        bm_txtable <- getBM(attributes,
                            filters="ensembl_transcript_id",
                            values=ensembl_transcript_ids,
                            mart=mart)
    }
    makeTranscriptDb(geneId = bm_txtable$ensembl_gene_id,
                     txId = bm_txtable$ensembl_transcript_id,
                     txChrom = bm_txtable$chromosome_name,
                     txStrand = ifelse(bm_txtable$strand == 1, "+", "-"),
                     txStart = bm_txtable$transcript_start,
                     txEnd = bm_txtable$transcript_end,
                     cdsStart = bm_txtable$cds_start,
                     cdsEnd = bm_txtable$cds_end,
                     exonStart = bm_txtable$exon_chrom_start,
                     exonEnd = bm_txtable$exon_chrom_end,
                     exonRank = bm_txtable$rank,
                     exonId = bm_txtable$ensembl_exon_id)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Comparing 2 TranscriptDb objects.
###

setMethod("as.data.frame", "TranscriptDb",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        COL2CLASS <- c(
            gene_id="character",
            tx_id="character",
            tx_name="character",
            tx_chrom="factor",
            tx_strand="factor",
            tx_start="integer",
            tx_end="integer",
            exon_rank="integer",
            exon_id="character",
            exon_chrom="factor",
            exon_strand="factor",
            exon_start="integer",
            exon_end="integer"
        )
        sql <- "SELECT
                  gene_id,
                  tx_id,
                  tx_name,
                  tx_chrom,
                  tx_strand,
                  tx_start, tx_end,
                  exon_rank,
                  exon_id,
                  exon_chrom,
                  exon_strand,
                  exon_start, exon_end,
                  CDS_id,
                  CDS_chrom,
                  CDS_strand,
                  CDS_start, CDS_end
                FROM transcripts
                  LEFT JOIN genes
                    ON (transcripts._tx_id=genes._tx_id)
                  INNER JOIN transcripts_rtree
                    ON (transcripts._tx_id=transcripts_rtree._tx_id)
                  INNER JOIN splicings
                    ON (transcripts._tx_id=splicings._tx_id)
                  INNER JOIN exons
                    ON (splicings._exon_id=exons._exon_id)
                  INNER JOIN exons_rtree
                    ON (exons._exon_id=exons_rtree._exon_id)
                  LEFT JOIN CDSs
                    ON (splicings._CDS_id=CDSs._CDS_id)
                  LEFT JOIN CDSs_rtree
                    ON (CDSs._CDS_id=CDSs_rtree._CDS_id)
                ORDER BY tx_chrom, tx_strand, tx_start, tx_end, tx_id, exon_rank"
        data <- dbGetQuery(x@conn, sql)
        setDataFrameColClass(data, COL2CLASS)
    }
)

compareTranscriptDbs <- function(txdb1, txdb2)
{
    if (!is(txdb1, "TranscriptDb")
     || !is(txdb2, "TranscriptDb"))
        stop("'txdb1' and 'txdb2' must be TranscriptDb objects")
    data1 <- as.data.frame(txdb1)
    data2 <- as.data.frame(txdb2)
    identical(data1, data2)
}

