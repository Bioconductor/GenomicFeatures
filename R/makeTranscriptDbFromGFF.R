### =========================================================================
### makeTranscriptDbFromGFF()
### -------------------------------------------------------------------------

## helper to clean up splicings
.filterDrop <- function(data, field){
  if(any(is.na(data[[field]]))){
    data <- data[!names(data) %in%  field]
  }
  data
}

.cleanSplicingsNAs <- function(data){
## check: exon_name, exon_chrom, exon_strand, cds_name
  fields <- c('exon_name','exon_chrom','exon_strand','cds_name')
  for(i in seq_len(length(fields))){
    if(any(colnames(data) %in% fields[i])){
      data <- .filterDrop(data, fields[i])
    }
  }
  data
}

## helper to clean up transcripts
.cleanTranscriptsNAs <- function(data){
## check: tx_name
  if(any(colnames(data) %in% 'tx_name')){
    data <- .filterDrop(data, 'tx_name')
  }
  data
}




## New strategy for getting the ranks: Use split on the starts and
## also on the strand to divide the pieces into chunks that can be
## computed on by mapply.  Then use mapply and call a function that
## considers both pieces together to make a decision and call order in
## the right way.  Then unlist and cbind etc.

## to be called by mapply
.assignRankings <- function(starts, strands){
  if(length(unique(strands)) >1 )
    stop("Exon rank inference cannot accomodate trans-splicing.")
  if (unique(strands) == "+") { 
    ord <- order(as.integer(starts))
  } else {
    ord <- order(as.integer(starts), decreasing=TRUE)
  }
  
  ## This new method will no longer sort the output (I don't expect it
  ## will matter.  But I will have to check...
 
  ord
}

.buildRanks <- function(exs){
    ## 1st we have to make 100% certain that we are we are sorted by tx_names
    exs <- exs[order(exs$tx_name),]
    ## then we can split by name
    starts <- split(exs[,"exon_start"], as.factor(exs$tx_name)[,drop=TRUE])
    strands <- split(exs[,"exon_strand"], as.factor(exs$tx_name)[,drop=TRUE])
    ranks <- unlist(mapply(.assignRankings, starts, strands))
    exs[,"exon_rank"] <- ranks
    exs
}



## Helper to deduce the rankings for each set of cds and exons...
.deduceExonRankings <- function(exs, format="gff"){
  message("Deducing exon rank from relative coordinates provided")
  ## And a warning for later (in case they were not watching)
  warning("Infering Exon Rankings.  If this is not what you expected, then please be sure that you have provided a valid attribute for exonRankAttributeName")
   res <- .buildRanks(exs)  
  ## then cast result to be data.frame 
  res <- data.frame(res, stringsAsFactors=FALSE)
  if(format=="gff"){
    colnames(res) <- c('exon_chrom','exon_start','exon_end','exon_strand',
                       'type','exon_name','tx_name','exon_rank')
  }else{
    colnames(res) <- c('tx_name','exon_rank','exon_chrom','exon_strand',
                       'exon_start','exon_end')
  }
  res$exon_start <- as.integer(res$exon_start)
  res$exon_end <- as.integer(res$exon_end)
  res$exon_rank <- as.integer(res$exon_rank)
  res
}


## Helpers to merge two frames together based only on ranges
## These two helpers stricly require columns:
## 'exon_chrom','exon_start','exon_end','exon_strand' OR
## 'cds_chrom','cds_start','cds_end','cds_strand' Along with 'tx_name' and
## 'exon_rank' (always) for each data.frame of input
.getUnusedColnames <- function(data){
  standardNames <- c('exon_chrom','exon_start','exon_end','exon_strand',
                     'cds_chrom','cds_start','cds_end','cds_strand')
  actualNames <- colnames(data)
  actualNames[!(actualNames %in% standardNames)]
}


## helper to massage 0 based rankings to be one based
.massageZeroBasedRankings <- function(splicings){
    if(any(splicings$exon_rank==0) && !any(splicings$exon_rank<0)){
        msg <- "Massaging the exon rank information so that it is one based counting instead of zero based."
        warning(paste(strwrap(msg, exdent=2), collapse="\n"),
                immediate.=TRUE, call.=FALSE)
        splicings$exon_rank <- splicings$exon_rank +1
    }
    splicings
}

## This helper might be generalizable for a common use case where people need
## to merge a pair of data.frames based on range data.  Though probably a
## different function that merges based on a pair of granges objects where you
## also want them to match based on a name (like a tx_name)
.mergeFramesViaRanges <- function(exs, cds){
  ##  we make GRanges objects so that we can range-match this stuff.
  exs <- exs[,!is.na(colnames(exs))]  
  cds <- cds[,!is.na(colnames(cds))]
  exsr <- GRanges(seqnames=Rle(exs$exon_chrom),
                  ranges=IRanges(start=exs$exon_start,end=exs$exon_end),
                  strand=exs$exon_strand,
                  exs[,.getUnusedColnames(exs)])
  cdsr <- GRanges(seqnames=Rle(cds$cds_chrom),
                  ranges=IRanges(start=cds$cds_start,end=cds$cds_end),
                  strand=cds$cds_strand,
                  cds[,.getUnusedColnames(cds)])
  ## call findOverlaps type='within' and have cdsr be the subject
  hits <- findOverlaps(query=cdsr,subject=exsr,type='within')
  if(any(duplicated(queryHits(hits)))){
    ## Get the names that go with each hit
    qTx <- mcols(cdsr[queryHits(hits)])$tx_name
    sTx <- mcols(exsr[subjectHits(hits)])$tx_name
    ## Then subset the hits to effectively filter based on names also matching
    hits <- hits[qTx == sTx]
  }
  
  ## To reassemble, I need to 1st get the matching bits:
  cdsExs <- cbind(exs[subjectHits(hits),],cds[queryHits(hits),])

  ## If there is info. for CDS rank and not for exons, then swap that
  ## over and send the user a message.
  if(all(is.na(exs$exon_rank)) && !all(is.na(cds$exon_rank))){
      msg <- "This file does not have information about exon rank, but it does have data for CDS rank.  Therefore we are applying the CDS rank information for the corresponding/overlapping exons.  This also means that any exons that did not have a CDS described will also not be present in the final TranscriptDb object."
      warning(paste(strwrap(msg, exdent=2), collapse="\n"),
              immediate.=TRUE, call.=FALSE)
      cdsExs$exon_rank <- cdsExs[[8]]
      splicings <- cdsExs
      splicings <- .massageZeroBasedRankings(splicings)
      return(splicings)
  }
  
  ## Finally I need to glue back the exon ranges that didn't have a cds...
  exsUnMatched <- exs[!(1:dim(exs)[1] %in% subjectHits(hits)),]
  emptys <- matrix(nrow = dim(exsUnMatched)[1],
                   ncol=dim(cds)[2]) ## ncol varies
  exsUn <- cbind(exsUnMatched, data.frame(emptys))
  names(exsUn) <- c(colnames(exs),colnames(cds))
  splicings <- rbind(cdsExs, exsUn)
  splicings <- .massageZeroBasedRankings(splicings)
  splicings ## coarse (contains everything that came in)
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the different data frames from GFF3
###   and return named list of tables
###

## helpers to pre-tidy the data
.checkExonRank <- function(data, gff, exonRankAttributeName){
  cnames <- colnames(data)
  if(!is.null(exonRankAttributeName)){
    exon_rank <- DataFrame(exonRankAttributeName=
                           mcols(gff)[[exonRankAttributeName]])
    data <- cbind(data,exon_rank)
    names(data) <- c(cnames, "exon_rank")
    
  }else{
    exon_rank <- DataFrame(rep(NA,length(start(gff))))
    data <- cbind(data,exon_rank)
    names(data) <- c(cnames, "exon_rank")    
  }
  data
}

.checkGeneIdAttrib <- function(data, gff, gffGeneIdAttributeName){
  if(!is.null(gffGeneIdAttributeName)){
    cnames <- colnames(data)
    gene_id <- DataFrame(gffGeneIdAttributeName=
                         mcols(gff)[[gffGeneIdAttributeName]])
    data <- cbind(data,gene_id)
    names(data) <- c(cnames, "gene_id")    
  }
  data
}

.prepareGFF3data.frame <- function(gff,exonRankAttributeName,
                                gffGeneIdAttributeName){
##   data <- DataFrame(seqnames=seqnames(gff),
##                     start=start(gff),
##                     end=end(gff),
##                     strand=strand(gff),
##                     type=mcols(gff)$type,
##                     ID=mcols(gff)$ID,
##                     Parent=mcols(gff)$Parent)
    
  dataR <- as(as(gff, "data.frame")[,1:5],"DataFrame")
  dataM <- as(mcols(gff), "DataFrame")
  if(dim(dataR)[1] == dim(dataM)[1]){
      data <- cbind(dataR,dataM) ## these will always be the same height
  }else{ stop("Error in .prepareGFF3data.frame, dataM and dataR should always match")}
  
  ## add ExonRank and geneID info if there is any
  data <- .checkExonRank(data, gff, exonRankAttributeName)
  data <- .checkGeneIdAttrib(data, gff, gffGeneIdAttributeName)
  
  ## Has a compressed col, so expand as needed.
  data <- expand(data, colnames="Parent", keepEmptyRows=TRUE )

  ## NOW convert data to a data.frame with no factors.
  data$seqnames <- as.character(data$seqnames)
  data$strand <- as.character(data$strand)
  data$type <- as.character(data$type)
  data$ID <- as.character(data$ID)
  data$Parent <- as.character(data$Parent)
  if("exonRankAttributeName" %in% colnames(data)){
      data$exonRankAttributeName <- as.integer(data$exonRankAttributeName)
  }
  if("gffGeneIdAttributeName" %in% colnames(data)){
      data$gffGeneIdAttributeName <- as.character(data$gffGeneIdAttributeName)
  }
  as.data.frame(data, stringsAsFactors=FALSE)
}


.prepareGFF3TXS <- function(data, useGenesAsTranscripts=FALSE){
  ## We absolutely require transcripts, genes and exons.
  message("extracting transcript information")
  if(useGenesAsTranscripts){
      txs <- data[data$type=="gene",]
  }else{
      txs <- data[data$type=="mRNA",]
  }
  if(dim(txs)[1] < 1){
      stop("No Transcript information found in gff file")
  }
  if(length(txs$ID) != length(unique(txs$ID))){
      stop("Unexpected transcript duplicates")}
  txs <- data.frame(txs, data.frame(tx_id=1:dim(txs)[1]),
                    stringsAsFactors=FALSE)
  as.data.frame(txs)
}


.prepareGFF3transcripts <- function(data,useGenesAsTranscripts){
  txs <- .prepareGFF3TXS(data,useGenesAsTranscripts)
  transcripts <- txs[,c("tx_id","ID","seqnames","strand","start","end")]
  names(transcripts) <- c("tx_id","tx_name","tx_chrom","tx_strand","tx_start",
                     "tx_end")
  ## Clean up any NA columns (any that are optional)
  transcripts <- .cleanTranscriptsNAs(transcripts)  
  as.data.frame(transcripts)
}

## This function probably has too many jobs.  Also it will need a refactor
.prepareGFF3genes <- function(data, transcripts, gffGeneIdAttributeName, gff,
                              useGenesAsTranscripts=FALSE){
  message("Extracting gene IDs")
  gns <- data[data$type=="gene",]
  if(dim(gns)[1] < 1){ ## that means there are no gene rows...
    ## Then we have to try and infer this from the transcript rows...
    if(!is.null(gffGeneIdAttributeName)){
      ## then try to compute gns using this other data...
      gns <- data.frame(tx_name=as.character(mcols(gff)$ID),
                        type=as.character(mcols(gff)$type),
                        gene_id=as.character(mcols(gff)[[gffGeneIdAttributeName]]),
                        stringsAsFactors=FALSE)
      gns <- gns[gns$type=="mRNA",]
      gns <- gns[,c("tx_name","gene_id")]
      ## Then merge in the more reliable tx_ids
      gns <- merge(gns, transcripts, by="tx_name")
      gns <- gns[,c("tx_id","gene_id")]
    }else{
      warning("No gene information present in gff file")
    }
  }else if(useGenesAsTranscripts){
      ## this case is mutually exclusive from the above case
      ## no need to warn twice
      txs <- suppressWarnings(.prepareGFF3TXS(data,
                                              useGenesAsTranscripts))
      gns <- txs[,c("tx_id","ID")] 
      names(gns) <- c("tx_id","gene_id") ## same as gns but with tx_id in tow.
  }else{
    if(length(gns$ID) != length(unique(gns$ID))){
      stop("Unexpected gene duplicates")}
    ## After testing for genes, here I get the actual data from mRNA rows...
    ## The only difference is that in this more normal case the Parents of
    ## these rows will be the genes that I detected previously.
    txs <- .prepareGFF3TXS(data,useGenesAsTranscripts)
    txsGene <- txs[,c("tx_id","Parent")] 
    names(txsGene) <- c("tx_id","gene_id")
    ## Then subset by gene_ids
    gns <- txsGene[txsGene$gene_id %in% gns$ID,]
  }
  as.data.frame(gns)
}


.prepareGFF3Fragments <- function(data, type){
    possibleCols <- c("seqnames","start","end","strand","type","ID","Parent",
                      "exon_rank","gene_id")
    expCols <- colnames(data) %in% possibleCols  
    res <- data[data$type==type,expCols]
    if (nrow(res) == 0L && type != "CDS")
        stop("No ", type, " information present in gff file")
    name <- paste0(tolower(type), "_id")
    colnames <- c("XXX_chrom","XXX_start","XXX_end","XXX_strand","type",
                  "XXX_name","tx_name","exon_rank")
    names(res) <- sub("XXX", tolower(type), colnames)
    res
}


.prepareGFF3Tables  <- function(gff,exonRankAttributeName,
                                gffGeneIdAttributeName,
                                useGenesAsTranscripts=FALSE){
  ## pre-clean the data
  data <- .prepareGFF3data.frame(gff,exonRankAttributeName,
                                gffGeneIdAttributeName)
  tables <- list()
  ## Get transcripts
  transcripts <- .prepareGFF3transcripts(data, useGenesAsTranscripts)
  tables[[1]] <- transcripts
  names(tables)[1] = "transcripts"
  ## Get genes
  tables[[2]] <- .prepareGFF3genes(data, transcripts,
                                   gffGeneIdAttributeName, gff,
                                   useGenesAsTranscripts)
  names(tables)[2] = "genes"

  message("Processing splicing information for gff3 file.")
  exs <- .prepareGFF3Fragments(data,type="exon")
  cds <- .prepareGFF3Fragments(data,type="CDS")
  
  ## if needed (usually needed for gff3), deduce the exon rankings from the
  ## order along the chromosome
  if(is.null(exonRankAttributeName)){
    exs <- .deduceExonRankings(exs, format="gff")
  }
  ## Then merge the two frames together based on range information
  splicings <- .mergeFramesViaRanges(exs, cds)
  
  ## now drop things and rename as needed
  splicings <- splicings[,c('exon_rank','exon_name','exon_chrom',
                            'exon_strand','exon_start','exon_end',
                            'cds_name','cds_start','cds_end','tx_name')]
  txIds <- unique(transcripts[,c("tx_id","tx_name")])
  splicings <- merge(txIds, splicings, by="tx_name")[,-1]
  ## Clean up any NA columns (any that are optional)
  splicings <- .cleanSplicingsNAs(splicings)
  
  if("cds_name" %in% colnames(splicings)){
  splicings$cds_name <- as(splicings$cds_name,"character")}
  if("cds_start" %in% colnames(splicings)){
      splicings$cds_start <- as(splicings$cds_start,"integer")}
  if("cds_end" %in% colnames(splicings)){
  splicings$cds_end <- as(splicings$cds_end,"integer")}
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

  message("Estimating transcript ranges.")
  ## for each transcript, we have to subset out the records and determine:
  ## start and stop based on strand using max and min.
  ## pre split the data for a substantial speedup
  sub <- as.data.frame(data, stringsAsFactors=FALSE)
  subs <- split(sub, as.factor(sub$transcript_id))
  ## and assign into a pre-allocated matrix
  res <- matrix(nrow = length(trns), ncol=6) ## always 6 wide

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


## helpers to pre-tidy the data
.prepareGTFdata.frame <- function(gff,exonRankAttributeName){
  data <- data.frame(seqnames=as.character(seqnames(gff)),
                       start=start(gff),
                       end=end(gff),
                       strand=as.character(strand(gff)),
                       type=as.character(mcols(gff)$type),
                       gene_id=as.character(mcols(gff)$gene_id),
                       transcript_id=as.character(mcols(gff)$transcript_id),
                       stringsAsFactors=FALSE)
  ## add ExonRank if there is any  
  if(!is.null(exonRankAttributeName)){
    data <- cbind(data,exon_rank=mcols(gff)[[exonRankAttributeName]])
  }else{
    data <- cbind(data,exon_rank=rep(NA,length(start(gff))))
  }
  data
}

## debug(GenomicFeatures:::.prepareGTFdata.frame)

.prepareGTFtranscripts <- function(data){
  ## We absolutely require transcripts, genes and exons.
  message("extracting transcript information")
  transcripts <- data
  if(length(unique(transcripts$transcript_id)) < 1){
    stop("No Transcript information present in gtf file")
  }else{
    ## GTF files require that we deduce the range of each transcript
    transcripts <- .deduceTranscriptsFromGTF(transcripts)
  }
  ## Clean up any NA columns (any that are optional)
  transcripts <- .cleanTranscriptsNAs(transcripts)
  transcripts
}

.prepareGTFgenes <- function(transcripts){
  message("Extracting gene IDs")
  gns <- transcripts[,c("tx_name","gene_id")]
  if(dim(gns)[1] < 1){warning("No gene information present in gtf file")
  }
  gns
}

.prepareGTFFragments <- function(data, type){
  res <- data[data$type==type,]
  res <- res[,c('transcript_id','exon_rank','seqnames','strand','start','end')]
  colnames <-  c('tx_name','exon_rank','XXX_chrom','XXX_strand',
                 'XXX_start','XXX_end')
  names(res) <- sub("XXX", tolower(type), colnames)
  res
}


## helper for preparing the GTF tables
.prepareGTFTables <- function(gff,exonRankAttributeName){
  ## pre-clean the data
  data <- .prepareGTFdata.frame(gff,exonRankAttributeName)
  tables <- list()
  ## get transcripts
  transcripts <- .prepareGTFtranscripts(data)
  ## get genes
  tables[[2]] <- .prepareGTFgenes(transcripts)
  names(tables)[2] <- "genes"

  ## once you have the genes, drop them from the transcripts
  transcripts <- transcripts[, c("tx_id","tx_name","tx_chrom","tx_strand",
                                 "tx_start","tx_end")]
  ## AND THEN assign them to the list
  tables[[1]] <- transcripts
  names(tables)[1] <- "transcripts"
  
  message("Processing splicing information for gtf file.")
  exs <- .prepareGTFFragments(data,type="exon")
  cds <- .prepareGTFFragments(data,type="CDS")
  
  ## if the exonRankAttributeName is not available ... then deduce.
  if(is.null(exonRankAttributeName)){
    exs <- .deduceExonRankings(exs,format="gtf")
  }
  
  ## no need to depend on having exon rank for cds too when we have this
  cdsExs <- .mergeFramesViaRanges(exs, cds)
  cdsExs <- cdsExs[,c("exon_rank","exon_chrom","exon_strand","exon_start",
                      "exon_end","cds_start","cds_end","tx_name")]
  
  ## now get the tx_ids by merging them in.
  txIds <- unique(transcripts[,c("tx_id","tx_name")])
  splicings <- merge(txIds, cdsExs, by="tx_name")[,-1]
  ## Clean up any NA columns (any that are optional)
  splicings <- .cleanSplicingsNAs(splicings)      
  tables[[3]] <- splicings
  names(tables)[3] <- "splicings"
  ## return all tables
  tables
}





## Helper to prepare the 'metadata' data frame.
## 

.prepareGFFMetadata <- function(file, dataSource=NA, species=NA,
                                      miRBaseBuild=NA)
{
    message("Prepare the 'metadata' data frame ... ",
            appendLF=FALSE)
    if (!isSingleStringOrNA(dataSource))
        stop("'dataSource' must be a a single string or NA")
    if (!isSingleStringOrNA(species))
        stop("'species' must be a a single string or NA")
    if (!isSingleStringOrNA(miRBaseBuild))
        stop("'miRBaseBuild' must be a a single string or NA")
    if (is.na(dataSource)) {
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
                   value=c(dataSource, species, miRBaseBuild)
                   )
    message("metadata: OK")
    metadata
}




### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### makeTranscriptDbFromGFF()
###

makeTranscriptDbFromGFF <- function(file,
                                    format=c("gff3", "gtf"),
                                    exonRankAttributeName=NULL,
                                    gffGeneIdAttributeName=NULL,
                                    chrominfo=NULL,
                                    dataSource=NA,
                                    species=NA,
                                    circ_seqs=DEFAULT_CIRC_SEQS,
                                    miRBaseBuild=NA,
                                    useGenesAsTranscripts=FALSE)
{
  format <- match.arg(format)
  
  ## start by importing the relevant features from the specified file
  feature.type <- c("gene", "mRNA", "exon", "CDS")
  gff <- import(file, format=format, feature.type=feature.type,
                asRangedData=FALSE)

  if(format=="gff3"){
    ## check that we have ID, Parent
    if(all(c("ID","Parent") %in% colnames(mcols(gff)))){
      tables <- .prepareGFF3Tables(gff, exonRankAttributeName,
                                   gffGeneIdAttributeName,
                                   useGenesAsTranscripts)
      ## results come back in list like: tables$transctripts etc.
    }
  }else if(format=="gtf"){
    ## check that we have gene_id and transcript_id
    if(all(c("gene_id","transcript_id")
           %in% colnames(mcols(gff)))){
      tables <- .prepareGTFTables(gff,exonRankAttributeName)
    }
  }
  ## TODO: verify that I have all I really need for metadata
  ## build up the metadata
  metadata <- .prepareGFFMetadata(file, dataSource, species, miRBaseBuild)

  ## If there is not chrominfo, then make one up best you can (no lengths)
  if(is.null(chrominfo)){
    message("Now generating chrominfo from available sequence names. No chromosome length information is available.")
    chroms <- unique(tables[["transcripts"]][["tx_chrom"]])
    chrominfo <- data.frame(chrom=chroms,
                            length=rep(NA,length(chroms)),
                            is_circular=matchCircularity(chroms, circ_seqs))
  }
  ## call makeTranscriptDb
  txdb <- makeTranscriptDb(transcripts=tables[["transcripts"]],
                           splicings=tables[["splicings"]],
                           genes=tables[["genes"]],
                           chrominfo=chrominfo,
                           metadata=metadata,
                           reassign.ids=TRUE)
  txdb
}

## ## TESTING GFF3
## gffFile=system.file("extdata","a.gff3")
## txdb <- makeTranscriptDbFromGFF(file=gffFile,
##                                 format="gff3",
##                                 dataSource="partial gtf file for Tomatoes
## donated anonymously for testing",
##                                 species="Solanum lycopersicum",
##                                 exonRankAttributeName = "nb_exon")
## saveDb(txdb,file="TESTGFF.sqlite")

## ## TESTING GTF
## gtfFile=system.file("extdata","Aedes_aegypti.partial.gtf")
## txdb <- makeTranscriptDbFromGFF(file=gtfFile,
##                                 format="gtf",
##                                 dataSource="ftp://ftp.ensemblgenomes.org/pub/metazoa/release-13/gtf/aedes_aegypti/",
##                                 species="Aedes aegypti")
## saveDb(txdb,file="TESTGTF.sqlite")


## TODO 5/3/12:
## )  Add some checks for columns in splicing and also transcripts for missing NA values, and drop if that is allowed
## ) Add unit tests
## ) fix any TODOs that still lie unanswered in this document.
## ) tidy the comments


##  library(GenomicFeatures);example(makeTranscriptDbFromGFF)

##  example(makeTranscriptDbFromGFF)


### Testing flybase file:
## flyFile = "dmel-4-r5.44.gff"
## txdb3 <- makeTranscriptDbFromGFF(file=flyFile,
##           format="gff3",
##           dataSource="gff file from flybase",
##           gffGeneIdAttributeName = "geneID",
##           species="Drosophila melanogaster")
## foo <- transcriptsBy(txdb3, by="gene")

## this file will not have any proper names etc.  Use it for testing.
## flyFile = "dmel-1000-r5.11.filtered.gff"


