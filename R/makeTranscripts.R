## Save method needs to connect to the DB, and write everything to a SQLite
## file It would be best if I could just take the existing handle, and then
## write that DB OUT.  I need to look a little closer to make this work.  Seth
## says there might be an easier way but that it would involve changes to
## RSQLite so perhaps I should hold off on this.

## Seth has also pointed out to me that I can pass a parameter into
## dbConnect() to set max.con = 100 so that the max number of connections is
## greater than 16).



saveFeatures <- function(file, annot){
  con <- annot@con
  ok <- sqliteCopyDatabase(con, file)
  stopifnot(ok)
}


getTranscriptData <- function(txAnnot){
  con <- txAnnot@con
  dbGetQuery(con, "SELECT * FROM transcripts")
}







.createAllDat <- function(frame){
  require(RSQLite)
  drv <- dbDriver("SQLite")  
  con <- dbConnect(drv)##, dbname="earlyTest.sqlite") ## **Temporarily** write to a file up front.
  sql <- "CREATE TABLE all_dat (
            _exon_id INTEGER PRIMARY KEY,
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
  dbGetQuery(con,sql)
  ##Just make sure that the order is correct...
  vals <- frame[,c("geneIds","ids","chrom","strand","txStart","txEnd","cdsStart",
                   "cdsEnd","exonStart","exonEnd","exonRank")]
  exon_ids <- 1:dim(frame)[1]
  vals <- cbind(exon_ids,vals)
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
  dbBeginTransaction(con)
  rset <- dbSendPreparedQuery(con, sqlIns, vals)
  
  dbClearResult(rset)
  dbCommit(con)
  
  dbGetQuery(con,"CREATE INDEX ad_tx_id on all_dat(tx_id)")
  dbGetQuery(con,"CREATE INDEX ad__exon_id on all_dat(_exon_id)")  
  con
}



.dropOldTable <- function(con, table){
  sql <- paste("DROP TABLE ",table,sep="")
  dbGetQuery(sql,con)
}


.createTranscriptsTable <- function(con){
  sql <- "CREATE TABLE transcripts (
            _tx_id INTEGER PRIMARY KEY,     --id a single transcript
            tx_id TEXT UNIQUE NOT NULL,     --text string for foreign transcript ID
            chromosome TEXT,                --redundant with info in exons
            strand TEXT )                   --redundant with info in exons
  "
  dbGetQuery(con,sql)

  sqli <- "INSERT INTO transcripts (tx_id, chromosome, strand) SELECT DISTINCT tx_id, chromosome, strand FROM all_dat"
  dbGetQuery(con,sqli)
  dbGetQuery(con,"CREATE INDEX t_tx_id on transcripts(tx_id)")
  dbGetQuery(con,"CREATE INDEX t__tx_id on transcripts(_tx_id)")
}


.createExonsTable <- function(con){
  sql <- "CREATE TABLE exons (
            _exon_id INTEGER PRIMARY KEY,   --id a single exon
            exon_id TEXT UNIQUE,   --text string for foreign exon ID
            chromosome TEXT,
            strand TEXT)
  "
  dbGetQuery(con,sql)

  sqli <- "INSERT INTO exons (chromosome, strand) SELECT chromosome, strand FROM all_dat"
  dbGetQuery(con,sqli)
  dbGetQuery(con,"CREATE INDEX e__exon_id on exons(_exon_id)")  
}


.createGenesTable <- function(con){
  sql <- "CREATE TABLE genes
    (gene_id VARCHAR(15),
     _tx_id VARCHAR(20),
     UNIQUE (gene_id,_tx_id),
     FOREIGN KEY (_tx_id) REFERENCES transcripts (_tx_id) )"
  dbGetQuery(con,sql)

  sqli <- "INSERT INTO genes SELECT DISTINCT ad.gene_id, t._tx_id
             FROM all_dat AS ad, transcripts AS t
             WHERE ad.tx_id = t.tx_id"
  dbGetQuery(con,sqli)
}


.createTranscriptTreeTable <- function(con){
  sql <- "CREATE VIRTUAL TABLE transcript_tree USING rtree
           (_tx_id INTEGER PRIMARY KEY,    --id a single transcript
            tx_start INTEGER,
            tx_end INTEGER,
            cds_start INTEGER,
            cds_end INTEGER)
            --FOREIGN KEY (_tx_id) REFERENCES  transcripts (_tx_id)"
  dbGetQuery(con,sql)

  sqli <- "INSERT INTO transcript_tree SELECT DISTINCT
             t._tx_id, ad.tx_start, ad.tx_end, ad.cds_start, ad.cds_end
             FROM transcripts AS t, all_dat AS ad
             WHERE ad.tx_id = t.tx_id"
  dbGetQuery(con,sqli)
}


##problem with new constraints means I have to make a "throw-away" table to ensure constraints!
.createExonTreePreTable <- function(con){
  sql <- "CREATE TABLE pre_exon_tree
           (_exon_id INTEGER PRIMARY KEY,    --id a single exon
            exon_start INTEGER,
            exon_end INTEGER,
            --UNIQUE (exon_start,exon_end),
            FOREIGN KEY (_exon_id) REFERENCES  exons (_exon_id) )"
  dbGetQuery(con,sql)

  sqli <- "INSERT INTO pre_exon_tree SELECT DISTINCT
             e._exon_id, ad.exon_start, ad.exon_end
             FROM exons AS e, all_dat AS ad
             WHERE ad._exon_id = e._exon_id"
  dbGetQuery(con,sqli)
}

.createExonTreeTable <- function(con){
  .createExonTreePreTable(con)
  sql <- "CREATE VIRTUAL TABLE exon_tree USING rtree
           (_exon_id INTEGER PRIMARY KEY,    --id a single exon
            exon_start INTEGER,
            exon_end INTEGER)
            --UNIQUE (exon_start,exon_end)
            --FOREIGN KEY (_exon_id) REFERENCES  exons (_exon_id)"
  dbGetQuery(con,sql)

  sqli <- "INSERT INTO exon_tree SELECT * FROM pre_exon_tree"
  dbGetQuery(con,sqli)
}





.createExonsTranscriptsTable <- function(con){
  sql <- "CREATE TABLE exons_transcripts (
            _exon_id INTEGER,            --id a single exon
            _tx_id INTEGER,              --id a single transcript
            exon_rank TEXT,
            UNIQUE(_exon_id,_tx_id),
            FOREIGN KEY (_exon_id) REFERENCES  exons  (_exon_id)
            FOREIGN KEY (_tx_id) REFERENCES  transcripts  (_tx_id))"
  dbGetQuery(con,sql)

  sqli <- "INSERT INTO exons_transcripts SELECT DISTINCT
             e._exon_id, t._tx_id, ad.exon_rank
             FROM exons as e, transcripts AS t, all_dat as ad
             WHERE ad.tx_id = t.tx_id AND e._exon_id = ad._exon_id"
  dbGetQuery(con,sqli)
}







## processFrame is an internal function to just get our code split up and
## formatted into our desired object.
.processFrame <- function(frame){

  ## Need to process this.  I want to jam it all into a DB
  con <- .createAllDat(frame)
  
  ##now just split things up to populate the other tables.
  .createTranscriptsTable(con)
  .createExonsTable(con)
  .createGenesTable(con)

  .createTranscriptTreeTable(con)

  .createExonsTranscriptsTable(con)
  .createExonTreeTable(con)


  ##drop the extra tables
##   .dropOldTable(con,"all_dat")
##   .dropOldTable(con,"pre_exon_tree")

  
  ## plan is to end with making the object
  new("TranscriptAnnotation",
      con = con)
}


















##TODO: These parameters need to be names that will identify the cols (not numbers).

convertExonsCommaSepFrame <- function(frame, exonColStart = "exonStarts", exonColEnd = "exonEnds"){

  eColStrtind = names(frame) %in% exonColStart
  eColEndind = names(frame) %in% exonColEnd
  
  ## Before I can call makeTranscripts, I have to split the exons up into their own rows.
  exonStart <- unlist(lapply(as.vector(frame[,eColStrtind]), function(x) strsplit(x,",")))
  exonEnd <- unlist(lapply(as.vector(frame[,eColEndind]), function(x) strsplit(x,",")))

  ## Need to also know how many times the rest needs to be expanded.
  lengths <- unlist(lapply(as.vector(frame[,eColStrtind]), function(x){ length(unlist(strsplit(x,","))) } ))

  ## I will use the subsetting shorthand to replicate the whole frame easily.
  ## This needs
  eColStrtEndind = !(names(frame) %in% c(exonColStart, exonColEnd))
  repCols = frame[rep(frame[,1],times=lengths ),eColStrtEndind]

  ## Derive the exonRanks:
  exonRank <- unlist(lapply(lengths, function(x) 1:x))
     
  data.frame(repCols, exonStart, exonEnd, exonRank, stringsAsFactors = FALSE)
  
}






makeTranscripts <- function(geneIds, ids, chrom, strand, txStart, txEnd, cdsStart,
                            cdsEnd, exonStart, exonEnd, exonRank,
                            exonId = vector(), exonColStart = "exonStarts",
                            exonColEnd = "exonEnds" ){
  if(length(exonId)<1){
    exonId = c(1:length(ids))
  }
  
  frame <- data.frame(geneIds = geneIds,
                      ids = ids,
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
  commas = FALSE
  if(length(grep(",", exonStart)) > 1 || length(grep(",", exonEnd)) > 1){
    commas = TRUE ## ie. we will be deriving the exonRanks
    frame = convertExonsCommaSepFrame(frame,
      exonColStart = exonColStart, exonColEnd = exonColEnd)
  }

  if(commas == FALSE){
    frame = cbind(frame, exonRank)
  }

  .processFrame(frame)
}


##TODO: The portion of this helper function that just retrieves org packages from
##common names could be generically useful in annotate.
.selectUCSCOrRefSeqMap <- function(type=c("knownGene","refGene"),
                       organism=c("human",
                                  "mouse",
                                  "rat",
                                  "dog",
                                  "chimp",
                                  "cow",
                                  "rhesus",
                                  "chicken",
                                  "zebrafish",
                                  "fly",
                                  "anopheles",
                                  "worm")){
  ##This method just applies some simple logic to decide which gene mapping
  ##to get and from which org packages.
  
  organism <- match.arg(organism)
  require(annotate)
  chip <- switch(organism,
                 human="org.Hs.eg",
                 mouse="org.Mm.eg",
                 rat="org.Rn.eg",
                 dog="org.Cf.eg",
                 chimp="org.Pt.eg",
                 cow="org.Bt.eg",
                 rhesus="org.Mmu.eg",
                 chicken="org.Gg.eg",
                 zebrafish="org.Dr.eg",
                 fly="org.Dm.eg",
                 anopheles="org.Ag.eg",
                 worm="org.Ce.eg")
  type <- match.arg(type)
  if(type=="knownGene"){
    map <- getAnnMap("UCSCKG", chip)
  }else if(type=="refGene"){
    map <- getAnnMap("REFSEQ", chip)
  }
  map
}


UCSCTranscripts <- function(type= "knownGene", genome="hg18",
                            organism="human",
                            exonColStart = "exonStarts",
                            exonColEnd = "exonEnds"){
  require("rtracklayer")
  session <- browserSession()
  genome(session) <- genome
  query <- ucscTableQuery(session, type)
  frame <- getTable(query)
  
  ## Mapping for UCSC refGene and knownGene types both use the same names for
  ## exonColStart and exonColEnd.  Params only exist as future insurance.

  ## For UCSC, we have to convert comma separated fields...
  frame <- convertExonsCommaSepFrame(frame,
    exonColStart = exonColStart, exonColEnd = exonColEnd)

  ##TODO: there are issues here with namespaces that keeps AnnotationDbi from
  ##being fully useful...

  ##Get the geneIDs and match them with the ids.  
  EGKGmap <- .selectUCSCOrRefSeqMap(type,organism)
  EGs <- unlist(mget(as.character(frame$name),revmap(EGKGmap), ifnotfound=NA)) 

  tx <- makeTranscripts(geneIds = EGs,
                        ids = frame$name,
                        chrom = frame$chrom,
                        strand = frame$strand,
                        txStart = frame$txStart,
                        txEnd = frame$txEnd,
                        cdsStart = frame$cdsStart, 
                        cdsEnd = frame$cdsEnd,
                        exonStart = frame$exonStart,
                        exonEnd = frame$exonEnd,
                        exonRank = frame$exonRank)  
}






## For people who want to tap biomaRt
## imagined usage:
## biomaRtTranscripts(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

BMTranscripts <- function(biomart="ensembl", dataset = "hsapiens_gene_ensembl"){
  require(biomaRt)
  mart = useMart(biomart=biomart, dataset=dataset)
  frame = getBM(mart=mart, attributes=c("ensembl_gene_id",
                             "ensembl_transcript_id",
                             "chromosome_name",
                             "strand",
                             "transcript_start",
                             "transcript_end",
                             "cds_start", ##TODO: cds_start frame of reference is WRONG! (++ to transcript_start)
                             "cds_end",
                             "ensembl_exon_id",
                             "exon_chrom_start",
                             "exon_chrom_end",
                             "rank" )) 
                             ##Add in missing values for cdsStart and cdsEnd AND gene_id

##   save(frame,file="BMFrame.Rda")
##   load("BMFrame.Rda")

  
  ##This one goes in the the next step the way that it came in.
  tx <- makeTranscripts(geneIds = frame$ensembl_gene_id,
                        ids = frame$ensembl_transcript_id,
                        chrom = frame$chromosome_name,
                        strand = frame$strand,
                        txStart = frame$transcript_start,
                        txEnd = frame$transcript_end,
                        cdsStart = frame$cds_start,
                        cdsEnd = frame$cds_end,
                        exonStart = frame$exon_chrom_start,
                        exonEnd = frame$exon_chrom_end,
                        exonRank = frame$rank,
                        exonId = frame$ensembl_exon_id) ##this is EXTRA info. that we don't have for UCSC - to make use of it we will want to attach it as an extra field in the DB.
 
}




##TODO: put a check in place to make double-damn sure that the unique tx_ids are equal in number to the number of unique exon IDs (and one per line) before attempting to make an all_dat table.  And if they are not, then thin them out...

##ALTERNATIVELY: I could just get the gene IDs separately for the BMFrame and add them on in R (instead of getting the whole thing from BM which results in the bloated frame)...



##TODO! Uncomment the constraint on exon starts and ends being the same, and change the way that these values are calculated so that they are NEVER assigned more than once.
