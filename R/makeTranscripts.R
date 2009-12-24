## Save method needs to connect to the DB, and write everything to a SQLite
## file It would be best if I could just take the existing handle, and then
## write that DB OUT.  I need to look a little closer to make this work.  Seth
## says there might be an easier way but that it would involve changes to
## RSQLite so perhaps I should hold off on this.

## Seth has also pointed out to me that I can pass a parameter into
## dbConnect() to set max.con = 100 so that the max number of connections is
## greater than 16).


### HP: I suggest to either rename saveFeatures()/loadFeatures() ->
### saveGenomicAnnotation()/loadGenomicAnnotation(), or to rename the
### GenomicAnnotation virtual class -> Features (or GenomicFeatures,
### probably better).
saveFeatures <- function(x, file) {
  conn <- x@conn
  ok <- sqliteCopyDatabase(conn, file)
  stopifnot(ok)
}


loadFeatures <- function(file) {
  if (!file.exists(file)) {
    stop("Cannot create a TranscriptAnnotation object without an actual database file.")
  }
  conn <- dbConnect(SQLite(), file)
  new("TranscriptAnnotation", conn=conn)
}


.createAllDat <- function(frame) {
  drv <- dbDriver("SQLite")  
  conn <- dbConnect(drv)##, dbname="earlyTest.sqlite") ## **Temporarily** write to a file up front.
  sql <- "CREATE TABLE all_dat (
            _exon_id INTEGER,
            gene_id VARCHAR(15),
            tx_id VARCHAR(15),
            chromosome VARCHAR(20),
            strand VARCHAR(3),
            tx_start VARCHAR(12),
            tx_end VARCHAR(12),
            cds_start VARCHAR(12),
            cds_end VARCHAR(12),
            exon_start VARCHAR(12),
            exon_end VARCHAR(12),
            exon_rank VARCHAR(12))"
  dbGetQuery(conn,sql)
  ##Just make sure that the order is correct...
  vals <- frame[,c("int_exon_id","geneId","txId","chrom",
                   "strand","txStart","txEnd","cdsStart","cdsEnd","exonStart",
                   "exonEnd","exonRank")]
##   exon_ids <- 1:dim(frame)[1]
##   vals <- cbind(exon_ids,vals)
  sqlIns <- "INSERT INTO all_dat (_exon_id,
                                  gene_id,
                                  tx_id,
                                  chromosome,
                                  strand,
                                  tx_start,
                                  tx_end,
                                  cds_start, 
                                  cds_end,
                                  exon_start,
                                  exon_end,
                                  exon_rank) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
  dbBeginTransaction(conn)
  rset <- dbSendPreparedQuery(conn, sqlIns, vals)
  
  dbClearResult(rset)
  dbCommit(conn)
  
  dbGetQuery(conn,"CREATE INDEX ad_tx_id on all_dat(tx_id)")
  dbGetQuery(conn,"CREATE INDEX ad__exon_id on all_dat(_exon_id)")  
  conn
}


.dropOldTable <- function(conn, table) {
  dbGetQuery(conn,paste("DROP TABLE ",table,sep=""))
}


.createTranscriptsTable <- function(conn) {
  sql <- "CREATE TABLE transcripts (
            _tx_id INTEGER PRIMARY KEY,     --id a single transcript
            tx_id TEXT UNIQUE NOT NULL,     --text string for foreign trnscpt ID
            chromosome TEXT,                --redundant with info in exons
            strand TEXT )                   --redundant with info in exons
  "
  dbGetQuery(conn,sql)

  sqli <- "INSERT INTO transcripts (tx_id, chromosome, strand) SELECT DISTINCT tx_id, chromosome, strand FROM all_dat"
  dbGetQuery(conn,sqli)
  dbGetQuery(conn,"CREATE INDEX t_tx_id on transcripts(tx_id)")
}


.createExonsTable <- function(conn) {
  sql <- "CREATE TABLE exons (
            _exon_id INTEGER PRIMARY KEY,   --id a single exon
            exon_id TEXT UNIQUE,   --text string for foreign exon ID
            chromosome TEXT,
            strand TEXT)
  "
  dbGetQuery(conn,sql)

  sqli <- "INSERT INTO exons (_exon_id, chromosome, strand) SELECT DISTINCT _exon_id, chromosome, strand FROM all_dat"
  dbGetQuery(conn,sqli)
  dbGetQuery(conn,"CREATE INDEX e__exon_id on exons(_exon_id)")  
}


.createGenesTable <- function(conn) {
  sql <- "CREATE TABLE genes
    (gene_id VARCHAR(15),
     _tx_id VARCHAR(20),
     UNIQUE (gene_id,_tx_id),
     FOREIGN KEY (_tx_id) REFERENCES transcripts (_tx_id) )"
  dbGetQuery(conn,sql)

  sqli <- "INSERT INTO genes SELECT DISTINCT ad.gene_id, t._tx_id
             FROM all_dat AS ad, transcripts AS t
             WHERE ad.tx_id = t.tx_id"
  dbGetQuery(conn,sqli)
}


.createTranscriptTreeTable <- function(conn) {
  sql <- "CREATE VIRTUAL TABLE transcript_tree USING rtree
           (_tx_id INTEGER PRIMARY KEY,    --id a single transcript
            tx_start INTEGER,
            tx_end INTEGER)
            --FOREIGN KEY (_tx_id) REFERENCES  transcripts (_tx_id)"
  dbGetQuery(conn,sql)

  sqli <- "INSERT INTO transcript_tree SELECT DISTINCT
             t._tx_id, ad.tx_start, ad.tx_end
             FROM transcripts AS t, all_dat AS ad
             WHERE ad.tx_id = t.tx_id"
  dbGetQuery(conn,sqli)
}


.createExonTreeTable <- function(conn) {
  sql <- "CREATE VIRTUAL TABLE exon_tree USING rtree
           (_exon_id INTEGER PRIMARY KEY,    --id a single exon
            exon_start INTEGER,
            exon_end INTEGER)
            --FOREIGN KEY (_exon_id) REFERENCES  exons (_exon_id)"
  dbGetQuery(conn,sql)

  sqli <- "INSERT INTO exon_tree SELECT DISTINCT
             e._exon_id, ad.exon_start, ad.exon_end
             FROM exons AS e, all_dat AS ad
             WHERE ad._exon_id = e._exon_id"
  dbGetQuery(conn,sqli)
}


.createExonsTranscriptsTable <- function(conn) {
  sql <- "CREATE TABLE exons_transcripts (
            _exon_id INTEGER,            --id a single exon
            _tx_id INTEGER,              --id a single transcript
            exon_rank INTEGER,
            UNIQUE(_exon_id,_tx_id),
            FOREIGN KEY (_exon_id) REFERENCES  exons  (_exon_id)
            FOREIGN KEY (_tx_id) REFERENCES  transcripts  (_tx_id))"
  dbGetQuery(conn,sql)

  sqli <- "INSERT INTO exons_transcripts SELECT DISTINCT
             e._exon_id, t._tx_id, ad.exon_rank
             FROM exons as e, transcripts AS t, all_dat as ad
             WHERE ad.tx_id = t.tx_id AND e._exon_id = ad._exon_id"
  dbGetQuery(conn,sqli)
}


## processFrame is an internal function to just get our code split up and
## formatted into our desired object.
.processFrame <- function(frame) {
  ## Need to process this.  I want to jam it all into a DB
  conn <- .createAllDat(frame)
  
  ##now just split things up to populate the other tables.
  .createTranscriptsTable(conn)
  .createExonsTable(conn)
  .createGenesTable(conn)
  .createTranscriptTreeTable(conn)
  .createExonsTranscriptsTable(conn)
  .createExonTreeTable(conn)

  ##drop the extra tables
  .dropOldTable(conn,"all_dat")
  
  ## plan is to end with making the object
  new("TranscriptAnnotation", conn=conn)
}


makeIdsForUniqueRows <- function(x, start="exonStart", end="exonEnd") {
    frame_names <- names(x)
    indices <- match(c("chrom", "strand", start, end), frame_names)
    if (any(is.na(indices)))
      stop("Some column names for x have not been specified correctly.")
    x <- x[,indices]
    ## NOTE: sorting is Locale specific, so different users will
    ## generate different IDs given the same 'x'.
    x_order <- do.call(order,x) 
    x_dups <- duplicated(x)
    ans <- vector("integer", length(x_order))
    ans[x_order] <- cumsum(!x_dups[x_order])
    ans
}


convertExonsCommaSepFrame <- function(frame, exonColStart = "exonStart",
                                      exonColEnd = "exonEnd")
{
    frame_names <- names(frame)
    startEndIndices <- match(c(exonColStart, exonColEnd), frame_names)
    if (any(is.na(startEndIndices)))
        stop("exonColStart or exonColEnd not found in column names")
    eColStrtind <- startEndIndices[1L]
    eColEndind <- startEndIndices[2L]
  
    ## Before I can call makeTranscripts, I have to split the exons up
    ## into their own rows.
    eStarts <- strsplit(as.character(frame[[eColStrtind]]), ",")
    lengths <- sapply(eStarts, length)
    exonStart <- unlist(eStarts)
    exonEnd <- unlist(strsplit(as.character(frame[[eColEndind]]), ","))

    ## I will use the subsetting shorthand to replicate the whole frame easily.
    ## This needs
    notStartEnd <- !(frame_names %in% c(exonColStart, exonColEnd))
    repCols <- frame[rep(seq_len(nrow(frame)), lengths), notStartEnd, drop = FALSE]
    keepNames <- names(repCols)
    ## Derive the exonRanks:
    exonRank <- unlist(lapply(lengths, function(x) seq_len(x)))

    ## FIXME: cleanup conversion of start/end to integer
    structure(data.frame(repCols, as.integer(exonStart), as.integer(exonEnd),
                         as.integer(exonRank), stringsAsFactors = FALSE,
                         row.names = NULL),
              names = c(keepNames, "exonStart", "exonEnd", "exonRank"))
}

### Creates a TranscriptAnnotation instance from vectors of data.
### All vectors must be of equal length.  The i-th element of the
### vectors represent data for a single exon. Transcript data will
### be repeated over all exons within the trancsript.
makeTranscripts <- function(geneId, txId, chrom, strand, txStart, txEnd,
                            cdsStart, cdsEnd, exonStart, exonEnd, exonRank,
                            exonId = vector()) {
  frame <- data.frame(geneId = geneId,
                      txId = txId,
                      chrom = chrom,
                      strand = strand,
                      txStart = txStart,
                      txEnd = txEnd,
                      cdsStart = cdsStart,
                      cdsEnd = cdsEnd,
                      exonStart = exonStart,
                      exonEnd = exonEnd,
                      stringsAsFactors = FALSE)
  
  ## Check if the exonEnd and exonStart are comma separated.
  if (length(grep(",", exonStart)) > 1 || length(grep(",", exonEnd)) > 1) {
    frame <- convertExonsCommaSepFrame(frame,
                                       exonColStart = "exonStart",
                                       exonColEnd = "exonEnd")
  } else {
    frame <- cbind(frame, as.integer(exonRank))
  }

  ##make unique IDs for unique exons.
  frame <- cbind(frame,
                 int_exon_id = makeIdsForUniqueRows(frame,
                                                    start = "exonStart",
                                                    end = "exonEnd"))
  .processFrame(frame)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### UCSCTranscripts() and UCSC-specific helper functions.
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

.makeTranscriptsFromUCSCTable <- function(txtable, track, organism)
{
    ## Map the transcript IDs in 'txtable$name' to the corresponding Entrez
    ## Gene IDs. Depending on the value of 'track' ("knownGene", "refGene"
    ## or "ensGene") this is done by using either the UCSCKG, the REFSEQ
    ## or the ENSEMBLTRANS map from the appropriate org package.
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
    map <- getAnnMap(mapname, orgpkg)
    ## Here we make the assumption that each transcript ID is mapped to at
    ## most one Entrez Gene ID, which will probably be true most of the
    ## time. For now we raise an error if this is not the case, but we might
    ## also want to handle this nicely in the future.
    tx2EG <- mget(as.character(txtable$name), revmap(map),
                  ifnotfound=NA_character_)
    if (any(sapply(tx2EG, length) > 1L))
        stop("some transcript IDs are mapped to multiple Entrez Gene IDs")
    EGs <- unname(unlist(tx2EG))  # can contain NAs
    txStart <- txtable$txStart + 1L
    cdsStart <- txtable$cdsStart + 1L
    exonStarts <- .shift.coordsInMultivaluedField(
                      as.character(txtable$exonStarts), 1L)
    makeTranscripts(geneId = EGs,
                    txId = txtable$name,
                    chrom = txtable$chrom,
                    strand = txtable$strand,
                    txStart = txStart,
                    txEnd = txtable$txEnd,
                    cdsStart = cdsStart,
                    cdsEnd = txtable$cdsEnd,
                    exonStart = exonStarts,
                    exonEnd = txtable$exonEnds)
}

UCSCTranscripts <- function(track=c("knownGene","refGene","ensGene"),
                            genome="hg18")
{
    track <- match.arg(track)
    if (!isSingleString(genome))
        stop("'genome' must be a single string")
    organism <- .getOrganismFromUCSCgenome(genome)
    if (track == "knownGene" && !(organism %in% c("Human", "Mouse", "Rat")))
        stop("UCSC \"knownGene\" track is only supported ",
             "for Human, Mouse and Rat")
    session <- browserSession()
    genome(session) <- genome
    query <- ucscTableQuery(session, track)
    table <- getTable(query)
    .makeTranscriptsFromUCSCTable(table, track, organism)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### BMTranscripts().
###
### For people who want to tap biomaRt.
### Typical usage:
###   BMTranscripts(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
###

BMTranscripts <- function(biomart="ensembl", dataset="hsapiens_gene_ensembl")
{
    mart <- useMart(biomart=biomart, dataset=dataset)
    txtable <- getBM(mart=mart,
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
    makeTranscripts(geneId = txtable$ensembl_gene_id,
                    txId = txtable$ensembl_transcript_id,
                    chrom = txtable$chromosome_name,
                    strand = txtable$strand,
                    txStart = txtable$transcript_start,
                    txEnd = txtable$transcript_end,
                    cdsStart = txtable$cds_start,
                    cdsEnd = txtable$cds_end,
                    exonStart = txtable$exon_chrom_start,
                    exonEnd = txtable$exon_chrom_end,
                    exonRank = txtable$rank,
                    exonId = txtable$ensembl_exon_id)
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
## tx = UCSCTranscripts()
## saveFeatures(tx, file="HG18.sqlite")
