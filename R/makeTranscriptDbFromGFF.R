### =========================================================================
### makeTranscriptDbFromGFF()
### -------------------------------------------------------------------------

## helper to calc the index
.computeStartInd <-function(i, res){
  ## The math on which indexes we need to fill with temp...
  if(i!=1){ ## normally we go to the 1st NA and add one
    startInd <- table(is.na(res[,1]))[["FALSE"]] + 1
  }else{startInd <- 1}
  startInd
}

## Helper to assign some rankings to the edge of a single set based on its
## strand information.
## This helper should look at the strand and then assign ranks one dir or the
## other.  If "-" then rank 1 is largest, if "+" then rank 1 is smallest.
.assignRankings <- function(dat){ # dat=es[[1]]
  dat <- as.matrix(dat)
  start <- "exon_start"
  rowLen <- dim(dat)[1]
  ## now sort (if needed, b/c I think normally this won't be necessary)
  if(dim(dat)[1]>1 && dat[,start][1] > dat[,start][rowLen]){
    ord <- order(dat[,start])
    dat <- dat[ord,]
  }
  ## now that dat will have been sorted (if it was needed), we can proceed to
  ## know that the 1st start value is the smallest number
  if(dat[,c("exon_strand")][[1]]=="+"){ #inference CANNOT HANDLE trans-splicing!
    dat[,c("exon_rank")] <- 1:rowLen
  }else{
    dat[,c("exon_rank")] <- rowLen:1
  }
  #as.matrix(dat) ## has to be a matrix to assign into one.
  dat
}


## Helper to deduce the rankings for each set of cds and exons...
.deduceExonRankings <- function(exs){
  message("Infering Exon Rankings.")
  res <- matrix(nrow = dim(exs)[1], ncol=9) ## all of it?
  ## split up the data
  es <- split(exs, as.factor(exs$tx_name))
  ## loop to assemble the result
  for(i in seq_len(length(es))){
    startInd <- .computeStartInd(i,res)
    endInd <- startInd + dim( es[[i]] )[1] - 1
    res[startInd:endInd,] <- .assignRankings(es[[i]])
  }
  ## then cast result to be data.frame 
  res <- data.frame(res, stringsAsFactors=FALSE)
  colnames(res) <- c('exon_chrom','exon_start','exon_end','exon_strand','type',
                     'exon_name','tx_name','exon_rank','exon_id')
  res$exon_start <- as.integer(res$exon_start)
  res$exon_end <- as.integer(res$exon_end)
  res$exon_rank <- as.integer(res$exon_rank)
  res$exon_id <- as.integer(res$exon_id)
  res
}




### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the different data frames from GFF3
###   and return named list of tables
### 

.prepareGFF3Tables <- function(gff,gffExonRankAttributeName){

  data <- DataFrame(seqnames=seqnames(gff),
                    start=start(gff),
                    end=end(gff),
                    strand=strand(gff),
                    type=gff$type,
                    ID=gff$ID,
                    Parent=gff$Parent,
                    exon_rank=gff[[gffExonRankAttributeName]])
  ## Assumption that we can rely on exon_rank being called nb_exon is
  ## probably FALSE.
  ## Therefore exon_rank will have to be generalized via argument.
  
  ## Has a compressed col, so expand as needed.
  data <- expand(data, colnames="Parent", keepEmptyRows=TRUE )

  ## NOW convert data to a data.frame with no factors.
  data$seqnames <- as.character(data$seqnames)
  data$strand <- as.character(data$strand)
  data$type <- as.character(data$type)
  data$ID <- as.character(data$ID)
  data$Parent <- as.character(data$Parent)
  data <- as.data.frame(data, stringsAsFactors=FALSE)
  
  tables <- list()
  ## We absolutely require transcripts, genes and exons.
  txs <- data[data$type=="mRNA",]
  if(dim(txs)[1] < 1){stop("No Transcript information present in gff file")
  }else{
    if(length(txs$ID) != length(unique(txs$ID))){
      stop("Unexpected transcript duplicates")}
    txs <- data.frame(txs, data.frame(tx_id=1:dim(txs)[1]),
                      stringsAsFactors=FALSE)
    txsSub <- txs[,c("tx_id","ID","seqnames","strand","start","end")]
    names(txsSub) <- c("tx_id","tx_name","tx_chrom","tx_strand","tx_start",
                       "tx_end")
    tables[[1]] <- as.data.frame(txsSub)
    names(tables)[1] = "transcripts"
  }

  gns <- data[data$type=="gene",]
  if(length(gns) < 1){stop("No gene information present in gff file")
  }else{
    if(length(gns$ID) != length(unique(gns$ID))){
      stop("Unexpected gene duplicates")}
    txsGene <- txs[,c("tx_id","Parent")]
    names(txsGene) <- c("tx_id","gene_id")
    ## Then subset by gene_ids
    gns <- txsGene[txsGene$gene_id %in% gns$ID,]
    tables[[2]] <- as.data.frame(gns)
    names(tables)[2] = "genes"
  }

  ## TODO, Too much repetition here: not mission critical but it bothers me to
  ## look at it.
  exs <- data[data$type=="exon",]
  if(length(exs) < 1){stop("No exon information present in gff file")
  }else{
    if(length(exs$ID) != length(unique(exs$ID))){
      stop("Unexpected exon duplicates")}
    exs <- data.frame(exs, data.frame(exon_id=1:dim(exs)[1]),
                      stringsAsFactors=FALSE)
    names(exs) <- c("exon_chrom","exon_start","exon_end","exon_strand","type",
                    "exon_name","tx_name","exon_rank","exon_id")
  }
  
  cds <- data[data$type=="CDS",]
  if(length(cds) < 1){warning("No CDS information present in gff file")
  }else{
    if(length(cds$ID) != length(unique(cds$ID))){
      stop("Unexpected cds duplicates")}
    cds <- data.frame(cds, data.frame(cds_id=1:dim(cds)[1]),
                      stringsAsFactors=FALSE)
    names(cds) <- c("cds_chrom","cds_start","cds_end","cds_strand","type",
                    "cds_name","tx_name","exon_rank","cds_id")    
  }
  

  ## For now lets AlWAYS deduce. TODO: make this optional in the event that
  ## exon_rank is provided by the file..
  ## deduce the exon rankings from the order along the chromosome
  exs <- .deduceExonRankings(exs)  
  ##  we have make GRanges objects so that we can range-match this stuff.
  exsr <- GRanges(seqnames=Rle(exs$exon_chrom),
                  ranges=IRanges(start=exs$exon_start,end=exs$exon_end),
                  strand=exs$exon_strand,
                  exs[,5:9])
  cdsr <- GRanges(seqnames=Rle(cds$cds_chrom),
                  ranges=IRanges(start=cds$cds_start,end=cds$cds_end),
                  strand=cds$cds_strand,
                  cds[,5:9])
  ## call findOverlaps type='within' and have cdsr be the subject
  hits <- findOverlaps(query=cdsr,subject=exsr,type='within')
  if(any(duplicated(queryHits(hits)))){
    ## Then we need to do filtering on all duplicated queryhits
    dupHits <- queryHits(hits[duplicated(queryHits(hits))])
    ## Get the names that go with each hit
    qTx <- values(cdsr[queryHits(hits)])$tx_name
    sTx <- values(exsr[subjectHits(hits)])$tx_name
    ## Then subset the hits
    hits <- hits[qTx == sTx]
  }
  
  ## To reassemble, I need to 1st get the matching bits:
  cdsExs <- cbind(exs[subjectHits(hits),],cds[queryHits(hits),])
  ## Now lets drop some of the cols
  cdsExs <- cdsExs[,c(1:9,11:12,15,18)]
  ## Finally I need to glue back the exon ranges that didn't have a cds...
  exsUnMatched <- exs[!(1:dim(exs)[1] %in% subjectHits(hits)),]
  emptys <- matrix(nrow = dim(exsUnMatched)[1], ncol=4)
  exsUn <- cbind(exsUnMatched, data.frame(emptys))
  names(exsUn) <- c(names(exsUnMatched),c('cds_start','cds_end','cds_name',
                                          'cds_id'))
  splicings <- rbind(cdsExs, exsUn)
  ## now drop things and rename as needed
  splicings <- splicings[,c('exon_rank','exon_id','exon_name','exon_chrom',
                            'exon_strand','exon_start','exon_end','cds_id',
                            'cds_name','cds_start','cds_end','tx_name')]
  txIds <- unique(txsSub[,c("tx_id","tx_name")])
  splicings <- merge(txIds, splicings, by="tx_name")[,-1]  
  tables[[3]] <- splicings
  names(tables)[3] <- "splicings"
  ## return all tables
  tables
}




### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the different data frames from GTF and return named list of tables 
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## helper for deducing the transcript range (max and min of all included exons)
.deduceTranscriptRangeData <- function(sub){
  ## Strand doesn't matter since the counting is always left to right.
  ## It's always the min of the starts and the max of the ends
  res <- c(unique(sub$transcript_id),
           unique(as.character(sub$seqnames)),
           unique(as.character(sub$strand)),
           min(sub$start),
           max(sub$end),
           unique(sub$gene_id) )
  res
}


## This works but is extremely slow.
.deduceTranscriptsFromGTF <- function(data){
  ## which transcripts?
  trns <- unique(data$transcript_id)

  message("Estimating transcript ranges - this might take a minute.")
  ## for each transcript, we have to subset out the records and determine:
  ## start and stop based on strand using max and min.
  res <- matrix(nrow = length(trns), ncol=6)
  ## pre split the data for a substantial speedup
  subs <- as.data.frame(data, stringsAsFactors=FALSE)
  subs <- split(subs, as.factor(subs$transcript_id))

  ## loop to assemble the result
  for(i in seq_len(length(trns))){
    res[i,] <- .deduceTranscriptRangeData(subs[[i]])
  }
  ## always make it a data.frame AFTER assiging in 
  res <- data.frame(res, stringsAsFactors=FALSE) 
  colnames(res) <- c("tx_name","tx_chrom","tx_strand","tx_start","tx_end",
                     "gene_id")
  ## start and end cannot be character() vectors
  res$tx_start <- as.integer(res$tx_start)
  res$tx_end <- as.integer(res$tx_end)
  
  ## Finally (at the end) assign the tx_ids to these.
  res <- data.frame(data.frame(tx_id=1:dim(res)[1]), res)
  res
}


## A faster version of exon/cds merging that takes the two data frames instead
## of two list objects
.fastMergeExonsCDSTables <- function(exs, cds){
  ## take advantage of the fact that we can just make a new ID col for each
  ## and merge on that...
  exsm <- data.frame(exs,
                mergeID=paste(exs$tx_name, exs$exon_rank, sep=""),
                stringsAsFactors=FALSE)
  cdsm <- data.frame(cds,
                mergeID=paste(cds$tx_name, cds$exon_rank, sep=""),
                stringsAsFactors=TRUE)[,c("cds_start","cds_end","mergeID")]
  
  res <- merge(exsm, cdsm, by="mergeID", all=TRUE)
  res <-res[,c("exon_rank","exon_chrom","exon_strand","exon_start",
               "exon_end","cds_start","cds_end","tx_name")]
  res  
}


.prepareGTFTables <- function(gff,gffExonRankAttributeName){

  ## Assemble the bits back together.
  data <- data.frame(seqnames=as.character(seqnames(gff)),
                     start=start(gff),
                     end=end(gff),
                     strand=as.character(strand(gff)),
                     type=as.character(gff$type),
                     gene_id=as.character(gff$gene_id),
                     transcript_id=as.character(gff$transcript_id),
                     exon_rank=gff[[gffExonRankAttributeName]],
                     stringsAsFactors=FALSE)
    
  tables <- list()
  ## We absolutely require transcripts, genes and exons.
  txs <- data  
  if(length(unique(txs$transcript_id)) < 1){
    stop("No Transcript information present in gtf file")
  }else{
    ## GTF files require that we deduce the range of each transcript
    system.time(txs <- .deduceTranscriptsFromGTF(txs))
  }

  gns <- txs[,c("tx_name","gene_id")]
  if(length(gns) < 1){stop("No gene information present in gtf file")
  }else{
    ## proceed to remove gene_ids column from txs (not needed there now)
    txs <- txs[, c("tx_id","tx_name","tx_chrom","tx_strand","tx_start",
                   "tx_end")]
  }
  tables[[1]] <- txs
  names(tables)[1] <- "transcripts"
  tables[[2]] <- gns
  names(tables)[2] <- "genes"

  message("Generating splicings from GTF file this will take some time.")
  exs <- data[data$type=="exon",]
  exs <- exs[,c('transcript_id','exon_rank','seqnames','strand','start','end')]
  names(exs) <- c('tx_name','exon_rank','exon_chrom','exon_strand',
                  'exon_start','exon_end')
  cds <- data[data$type=="CDS",]
  cds <- cds[,c('transcript_id','exon_rank','start','end')]
  names(cds) <- c('tx_name','exon_rank','cds_start','cds_end')

  ## make pre-split versions of each table
  es <- split(exs, as.factor(exs$tx_name))
  cs <- split(cds, as.factor(cds$tx_name))
  ## cdsExs <- .mergeExonsCDSTables(es, cs)
  cdsExs <- .fastMergeExonsCDSTables(exs, cds)
  
  ## now get the tx_ids by merging them in.
  txIds <- unique(txs[,c("tx_id","tx_name")])
  splicings <- merge(txIds, cdsExs, by="tx_name")[,-1]
  tables[[3]] <- splicings
  names(tables)[3] <- "splicings"
  ## return all tables
  tables
}


## Helper to prepare the 'metadata' data frame.
## 

.prepareGFFMetadata <- function(dataSource,species,miRBaseBuild)
{
    message("Prepare the 'metadata' data frame ... ",
            appendLF=FALSE)
    if(is.null(miRBaseBuild)){ miRBaseBuild <- NA }
    if(missing(dataSource)){ miRBaseBuild <- "GFF or GTF file" }
    metadata <- data.frame(
                   name=c("Data source",
                     "Genus and Species",
                     "miRBase build ID"),
                   value=c(dataSource,
                     species,
                     miRBaseBuild)
                   )
    message("metadata: OK")
    metadata
}





### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTranscriptDbFromGFF()
###

makeTranscriptDbFromGFF <- function(file,
                                    format=c("gff3", "gtf"),
                                    gffExonRankAttributeName,
                                    dataSource,
                                    species,
                                    circ_seqs=DEFAULT_CIRC_SEQS,
                                    miRBaseBuild=NULL)
{
  format <- match.arg(format)
  ## start by importing the file
  gff <- import(file, format=format)

  if(format=="gff3"){
    ## check that we have ID, Parent and nb_exon???
    if(all(c("ID","Parent", gffExonRankAttributeName) %in% colnames(gff))){
      tables <- .prepareGFF3Tables(gff, gffExonRankAttributeName)
      ## results come back in list like: tables$transctripts etc.
    }
  }else if(format=="gtf"){
    if(missing(gffExonRankAttributeName)){
      gffExonRankAttributeName  <- "exon_number" }
    ## check that we have gene_id and transcript_id
    if(all(c("gene_id","transcript_id",gffExonRankAttributeName)
           %in% colnames(gff))){
      tables <- .prepareGTFTables(gff,gffExonRankAttributeName)
    }
  }
  ## TODO: add arguments to allow construction of rudimentary metadata
  metadata <- .prepareGFFMetadata(dataSource,species,miRBaseBuild)

  ## BUT if there isn't one, lets make an NA filled chrominfo like this:
  chroms <- unique(tables[["transcripts"]][["tx_chrom"]])
  chrominfo <- data.frame(chrom=chroms, length=rep(NA,length(chroms)),
             is_circular=GenomicFeatures:::matchCircularity(chroms, circ_seqs))
  ## call makeTranscriptDb
  txdb <- makeTranscriptDb(transcripts=tables[["transcripts"]],
                           splicings=tables[["splicings"]],
                           genes=tables[["genes"]],
                           chrominfo=chrominfo,
                           metadata=metadata)
  txdb
}









## ## TESTING GFF3
## gffFile=system.file("extdata","a.gff3")
## txdb <- makeTranscriptDbFromGFF(file=gffFile,
##                                 format="gff3",
##                                 dataSource="partial gtf file for Tomatoes
## donated anonymously for testing",
##                                 species="Solanum lycopersicum",
##                                 gffExonRankAttributeName = "nb_exon",
##                                 circ_seqs=DEFAULT_CIRC_SEQS,
##                                 miRBaseBuild=NULL)
## saveDb(txdb,file="TESTGFF.sqlite")



## ## TESTING GTF
## gtfFile=system.file("extdata","Aedes_aegypti.partial.gtf")
## source("../makeTranscriptDbFromGTF.R")
## txdb <- makeTranscriptDbFromGFF(file=gtfFile,
##                                 format="gtf",
##                                 dataSource="ftp://ftp.ensemblgenomes.org/pub/metazoa/release-13/gtf/aedes_aegypti/",
##                                 species="Aedes aegypti",
##                                 circ_seqs=DEFAULT_CIRC_SEQS,
##                                 miRBaseBuild=NULL)
## saveDb(txdb,file="TESTGTF.sqlite")







## TODO 5/3/12:
## ) alter code for gff parsing so that it can work if there is an exon rank supplied and add code to gtf parsing so that it can infer the ranks.  (right now both of these are separated.  Basically, generalize the range matching strategy and use it whenever inference is required, meanwhile whenever ranges are provided, they might only exist for exons, so you should also use the range match strategy to finish in that case as well.  So the question is just one of whether or not we have exons and whether or not we have to call the helpers to infer them (and also standardizing our column names earlier).  This refactor will also reduce the amount of code in this document.
## ) Add unit tests
## ) fix any TODOs that still lie unanswered in this document.
