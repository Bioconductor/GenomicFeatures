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
    conn <- x@conn
    ok <- sqliteCopyDatabase(conn, file)
    if (!ok)
        stop("failed to write 'x' to file '", file, "'")
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
### makeTranscriptDb().
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
                                    chrom,
                                    strand,
                                    start,
                                    end,
                                    colnames)
{
    table <- data.frame(
        internal_id=internal_id,
        external_id=external_id,
        chrom=chrom,
        strand=strand,
        start=start,
        end=end,
        stringsAsFactors=FALSE)
    table <- unique(table)
    ## Create the '<feature>s' table.
    if (all(is.na(external_id)))
        null_constraint <- ""
    else
        null_constraint <- " NOT NULL"
    sql <- c(
        "CREATE TABLE ", feature, "s (\n",
        "  ", colnames[1L], " INTEGER PRIMARY KEY,\n",
        "  ", colnames[2L], " TEXT UNIQUE", null_constraint, ",\n",
        "  ", colnames[3L], " TEXT NOT NULL,\n",
        "  ", colnames[4L], " TEXT NOT NULL\n",
        ")")
    res <- dbSendQuery(conn, paste(sql, collapse=""))
    dbClearResult(res)
    ## Fill the '<feature>s' table.
    ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
    ## when the nb of rows to insert is 0, hence the following test.
    if (nrow(table) != 0L) {
        sql <- c("INSERT INTO ", feature, "s VALUES (?,?,?,?)")
        dbBeginTransaction(conn)
        res <- dbSendPreparedQuery(conn, paste(sql, collapse=""),
                                   table[1:4])
        dbClearResult(res)
        dbCommit(conn)
    }
    ## Create the '<feature>s_rtree' table.
    sql <- c(
        "CREATE VIRTUAL TABLE ", feature, "s_rtree USING rtree (\n",
        "  ", colnames[1L], " INTEGER PRIMARY KEY,\n",
        "  ", colnames[5L], " INTEGER, -- NOT NULL is implicit in rtree\n",
        "  ", colnames[6L], " INTEGER  -- NOT NULL is implicit in rtree\n",
        "  -- FOREIGN KEY (", colnames[1L], ") REFERENCES ", feature, "s\n",
        ")")
    res <- dbSendQuery(conn, paste(sql, collapse=""))
    dbClearResult(res)
    ## Fill the '<feature>s_rtree' table.
    if (nrow(table) != 0L) {
        sql <- c("INSERT INTO ", feature, "s_rtree VALUES (?,?,?)")
        dbBeginTransaction(conn)
        res <- dbSendPreparedQuery(conn, paste(sql, collapse=""),
                                   table[c(1L, 5:6)])
        dbClearResult(res)
        dbCommit(conn)
    }
}

### TODO: Call the table 'transcripts_exons' and define its cols in the
### following order: '_tx_id', 'exon_rank', '_exon_id'.
.writeExonsTranscriptsTable <- function(conn,
                                        internal_exon_id,
                                        internal_tx_id,
                                        exonRank)
{
    table <- data.frame(
        internal_exon_id=internal_exon_id,
        internal_tx_id=internal_tx_id,
        exonRank=exonRank,
        stringsAsFactors=FALSE)
    table <- unique(table)
    ## Create the 'exons_transcripts' table.
    sql <- c(
        "CREATE TABLE exons_transcripts (\n",
        "  _exon_id INTEGER NOT NULL,\n",
        "  _tx_id INTEGER NOT NULL,\n",
        "  exon_rank INTEGER NOT NULL,\n",
        "  UNIQUE (_tx_id, exon_rank),\n",
        "  FOREIGN KEY (_exon_id) REFERENCES exons,\n",
        "  FOREIGN KEY (_tx_id) REFERENCES transcripts\n",
        ")")
    res <- dbSendQuery(conn, paste(sql, collapse=""))
    dbClearResult(res)
    ## Fill the 'exons_transcripts' table.
    ## sqliteExecStatement() (SQLite backend for dbSendPreparedQuery()) fails
    ## when the nb of rows to insert is 0, hence the following test.
    if (nrow(table) != 0L) {
        sql <- "INSERT INTO exons_transcripts VALUES (?,?,?)"
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

### Creates a TranscriptDb instance from vectors of data.
### All vectors must be of equal length.  The i-th element of the
### vectors represent data for a single exon. Transcript data will
### be repeated over all exons within the trancsript.
makeTranscriptDb <- function(geneId, txId, chrom, strand, txStart, txEnd,
                             cdsStart, cdsEnd, exonStart, exonEnd, exonRank,
                             exonId=NULL)
{
    if (!is.character(geneId))
        stop("'geneId' must be a character vector")
    if (!is.character(txId) || any(is.na(txId)))
        stop("'txId' must be a character vector with no NAs")
    chrom <- .argAsCharacterFactorWithNoNAs(chrom, "chrom")
    strand <- .argAsCharacterFactorWithNoNAs(strand, "strand")
    txStart <- .argAsIntegerWithNoNAs(txStart, "txStart")
    txEnd <- .argAsIntegerWithNoNAs(txEnd, "txEnd")
    cdsStart <- .argAsIntegerWithNoNAs(cdsStart, "cdsStart")
    cdsEnd <- .argAsIntegerWithNoNAs(cdsEnd, "cdsEnd")
    exonStart <- .argAsIntegerWithNoNAs(exonStart, "exonStart")
    exonEnd <- .argAsIntegerWithNoNAs(exonEnd, "exonEnd")
    exonRank <- .argAsIntegerWithNoNAs(exonRank, "exonRank")
    if (is.null(exonId)) {
        exonId <- rep.int(as.character(NA), length(exonStart))
    } else {
        if (!is.character(exonId) || any(is.na(exonId)))
            stop("when supplied, 'exonId' must be a character vector ",
                 "with no NAs")
    }
    internal_tx_id <- .makeInternalIdsFromExternalIds(txId) 
    internal_exon_id <- .makeInternalIdsForUniqueLocs(chrom, strand,
                                                      exonStart, exonEnd)

    conn <- dbConnect(SQLite(), dbname="") ## we'll write the db to a temp file
    .writeFeatureCoreTables(conn, "transcript",
        internal_tx_id, txId, chrom, strand, txStart, txEnd,
        c("_tx_id", "tx_id", "chromosome", "strand",
          "tx_start", "tx_end"))
    .writeFeatureCoreTables(conn, "exon",
        internal_exon_id, exonId, chrom, strand, exonStart, exonEnd,
        c("_exon_id", "exon_id", "chromosome", "strand",
          "exon_start", "exon_end"))
    .writeExonsTranscriptsTable(conn,
        internal_exon_id, internal_tx_id, exonRank)
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
    txtable <- data.frame(
        geneId=rep.int(EGs, exon_count),
        txId=rep.int(ucsc_txtable$name, exon_count),
        chrom=rep.int(ucsc_txtable$chrom, exon_count),
        strand=rep.int(ucsc_txtable$strand, exon_count),
        txStart=rep.int(ucsc_txtable$txStart, exon_count) + 1L,
        txEnd=rep.int(ucsc_txtable$txEnd, exon_count),
        cdsStart=rep.int(ucsc_txtable$cdsStart, exon_count) + 1L,
        cdsEnd=rep.int(ucsc_txtable$cdsEnd, exon_count),
        exonStart=unlist(exonStarts) + 1L,
        exonEnd=unlist(exonEnds),
        exonRank=IRanges:::mseq(rep.int(1L, length(exon_count)), exon_count),
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
###

makeTranscriptDbFromBiomart <- function(biomart="ensembl",
                                        dataset="hsapiens_gene_ensembl")
{
    mart <- useMart(biomart=biomart, dataset=dataset)
    bm_txtable <- getBM(mart=mart,
                     attributes=c("ensembl_gene_id",
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
                                  "rank"))
    makeTranscriptDb(geneId = bm_txtable$ensembl_gene_id,
                     txId = bm_txtable$ensembl_transcript_id,
                     chrom = bm_txtable$chromosome_name,
                     strand = bm_txtable$strand,
                     txStart = bm_txtable$transcript_start,
                     txEnd = bm_txtable$transcript_end,
                     cdsStart = bm_txtable$cds_start,
                     cdsEnd = bm_txtable$cds_end,
                     exonStart = bm_txtable$exon_chrom_start,
                     exonEnd = bm_txtable$exon_chrom_end,
                     exonRank = bm_txtable$rank,
                     exonId = bm_txtable$ensembl_exon_id)
    ##this is EXTRA info. that we don't have for UCSC - to make use of it we
    ##will want to attach it as an extra field in the DB.
}


##TODO: put a check in place to make double-damn sure that the unique tx_ids
##are equal in number to the number of unique exon IDs (and one per line)
##before attempting to make an all_dat table.  And if they are not, then thin
##them out...

##ALTERNATIVELY: I could just get the gene IDs separately for the BMFrame and
##add them on in R (instead of getting the whole thing from BM which results
##in the bloated frame)...


##The trouble with the BMFrame Data is that the data have variable cds start
##and ends for a given tx-id.  So when I go to insert this into the DB, it
##violates a constraint in place to insure row-wise integrity for each tx
##combined with the constraint on only allowing one tx_id per row (table
##constraint).

##I have already considered ordering based on the cds start and ends when
##assigning tx_ids.  But this is probably a bad idea because we don't want to
##really even get into that whole game and we kind of need to trust the
##incoming data a certain amount.

##ALSO, it appears that the cds_start and cds_end data from biomart are not
##only counted based on a transcript frame of reference, but that they are
##also split up based on exons...  Which complicates things in terms of
##pre-processing.


##And actually, the cds_start and cds_end from UCSC are also kind of screwy
##(they simply look WRONG).


##TODO: These frames all need to be adjusted so that their counts are correct
##(pre-processed to adjust to the same counting offset as we use in R).


##used to make the testDB:
##Insert the first bit into UCSCTranscript
##frame = frame[c(1:500, 100000:100100, 120000:120100, 150000:150100,
##170000:170100, 200000:200100),]
## library(GenomicFeatures)
## txdb = makeTranscriptDbFromUCSC()
## saveFeatures(txdb, file="HG18.sqlite")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Comparing 2 TranscriptDb objects.
###

setMethod("as.data.frame", "TranscriptDb",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        COL2CLASS <- c(
            gene_id="character",
            tx_id="character",
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
                  transcripts.chromosome AS tx_chrom,
                  transcripts.strand AS tx_strand,
                  tx_start, tx_end,
                  exon_rank,
                  exon_id,
                  exons.chromosome AS exon_chrom,
                  exons.strand AS exon_strand,
                  exon_start, exon_end
                FROM transcripts
                  LEFT JOIN genes
                    ON (transcripts._tx_id=genes._tx_id)
                  INNER JOIN transcripts_rtree
                    ON (transcripts._tx_id=transcripts_rtree._tx_id)
                  INNER JOIN exons_transcripts
                    ON (transcripts._tx_id=exons_transcripts._tx_id)
                  INNER JOIN exons
                    ON (exons_transcripts._exon_id=exons._exon_id)
                  INNER JOIN exons_rtree
                    ON (exons._exon_id=exons_rtree._exon_id)
                ORDER BY tx_id, exon_rank"
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

