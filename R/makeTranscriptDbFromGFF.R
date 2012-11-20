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


## helper to calc the index
.computeStartInd <-function(i, res){
  ## The math on which indexes we need to fill with temp...
  if(i!=1){ ## normally we go to the 1st NA and add one
    startInd <- table(is.na(res[,1]))[["FALSE"]] + 1
  }else{startInd <- 1}
  startInd
}

## Helper to assign some rankings to the edge of a single set based on the
## strand information.
## This helper should look at the strand and then assign ranks one dir or the
## other.  If "-" then rank 1 is largest, if "+" then rank 1 is smallest.
.assignRankings <- function(dat){ # dat=es[[1]]
  dat <- as.matrix(dat)
  start <- dat[,"exon_start"]
  strand <- dat[,"exon_strand"]
  if(length(unique(strand)) >1 )
    stop("Inference CANNOT HANDLE trans-splicing.")
  if (strand[1]=="+") {
    ord <- order(as.integer(start))
  } else {
    ord <- order(as.integer(start), decreasing=TRUE)
  }
  ## now sort (needed or not (cheap enough!))
  dat <- dat[ord,,drop=FALSE]
  dat[,"exon_rank"] <- seq_len(nrow(dat))
  dat
}


.buildRanks <- function(es, res){
  ## loop to assemble the result
  for(i in seq_len(length(es))){
    startInd <- .computeStartInd(i,res)
    endInd <- startInd + dim( es[[i]] )[1] - 1
    res[startInd:endInd,] <- .assignRankings(es[[i]])
  }
  res
}

## Helper to deduce the rankings for each set of cds and exons...
.deduceExonRankings <- function(exs, format="gff"){
  message("Deducing exon rank from relative coordinates provided")
  ## And a warning for later (in case they were not watching)
  warning("Infering Exon Rankings.  If this is not what you expected, then please be sure that you have provided a valid attribute for exonRankAttributeName")
  res <- matrix(nrow = dim(exs)[1], ncol=dim(exs)[2]) ## ncol=9?  
  ## split up the data
  es <- split(exs, as.factor(exs$tx_name)[,drop=TRUE])
  ## loop to assemble the result
  res <- .buildRanks(es, res)  
  ## then cast result to be data.frame 
  res <- data.frame(res, stringsAsFactors=FALSE)
  if(format=="gff"){
    colnames(res) <- c('exon_chrom','exon_start','exon_end','exon_strand',
                       'type','exon_name','tx_name','exon_rank','exon_id')
    res$exon_id <- as.integer(res$exon_id)
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

## This helper might be generalizable for a common use case where people need
## to merge a pair of data.frames based on range data.  Though probably a
## different function that merges based on a pair of granges objects where you
## also want them to match based on a name (like a tx_name)
.mergeFramesViaRanges <- function(exs, cds){
  ##  we have make GRanges objects so that we can range-match this stuff.
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
  
  ## Finally I need to glue back the exon ranges that didn't have a cds...
  exsUnMatched <- exs[!(1:dim(exs)[1] %in% subjectHits(hits)),]
  emptys <- matrix(nrow = dim(exsUnMatched)[1],
                   ncol=dim(cds)[2]) ## ncol varies
  exsUn <- cbind(exsUnMatched, data.frame(emptys))
  names(exsUn) <- c(colnames(exs),colnames(cds))
  splicings <- rbind(cdsExs, exsUn)
  splicings ## coarse (contains everything that came in)
}



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Extract the different data frames from GFF3
###   and return named list of tables
###

## helpers to pre-tidy the data
.checkExonRank <- function(data, gff, exonRankAttributeName){
  if(!is.null(exonRankAttributeName)){
    exon_rank <- DataFrame(gff[[exonRankAttributeName]])
    data <- cbind(data,exon_rank)
  }else{
    exon_rank <- DataFrame(rep(NA,length(start(gff))))
    data <- cbind(data,exon_rank)
  }
  data
}

.checkGeneIdAttrib <- function(data, gff, gffGeneIdAttributeName){
  if(!is.null(gffGeneIdAttributeName)){
    gene_id <- DataFrame(gff[[gffGeneIdAttributeName]])
    data <- cbind(data,gene_id)
  }
  data
}

.prepareGFF3data.frame <- function(gff,exonRankAttributeName,
                                gffGeneIdAttributeName){
  data <- DataFrame(seqnames=seqnames(gff),
                    start=start(gff),
                    end=end(gff),
                    strand=strand(gff),
                    type=gff$type,
                    ID=gff$ID,
                    Parent=gff$Parent)
  ## add ExonRank and geneID info if there is any
  data <- .checkExonRank(data, gff, exonRankAttributeName)
  data <- .checkGeneIdAttrib(data, gff, exonRankAttributeName)
  
  ## Has a compressed col, so expand as needed.
  data <- expand(data, colnames="Parent", keepEmptyRows=TRUE )

  ## NOW convert data to a data.frame with no factors.
  data$seqnames <- as.character(data$seqnames)
  data$strand <- as.character(data$strand)
  data$type <- as.character(data$type)
  data$ID <- as.character(data$ID)
  data$Parent <- as.character(data$Parent)
  as.data.frame(data, stringsAsFactors=FALSE)
}


.prepareGFF3TXS <- function(data){
  ## We absolutely require transcripts, genes and exons.
  message("extracting transcript information")
  txs <- data[data$type=="mRNA",]
  if(dim(txs)[1] < 1){stop("No Transcript information present in gff file")
  }else{
    if(length(txs$ID) != length(unique(txs$ID))){
      stop("Unexpected transcript duplicates")}
    txs <- data.frame(txs, data.frame(tx_id=1:dim(txs)[1]),
                      stringsAsFactors=FALSE)
  }
  as.data.frame(txs)
}


.prepareGFF3transcripts <- function(data){
  txs <- .prepareGFF3TXS(data)
  transcripts <- txs[,c("tx_id","ID","seqnames","strand","start","end")]
  names(transcripts) <- c("tx_id","tx_name","tx_chrom","tx_strand","tx_start",
                     "tx_end")
  ## Clean up any NA columns (any that are optional)
  transcripts <- .cleanTranscriptsNAs(transcripts)  
  as.data.frame(transcripts)
}

.prepareGFF3genes <- function(data, transcripts, gffGeneIdAttributeName, gff){
  message("Extracting gene IDs")
  gns <- data[data$type=="gene",]
  if(dim(gns)[1] < 1){
    ## Then we have to try and infer this from the transcript rows...
    if(!is.null(gffGeneIdAttributeName)){
      ## then try to compute gns using this other data...
      gns <- data.frame(tx_name=as.character(gff$ID),
                        type=as.character(gff$type),
                        gene_id=as.character(gff[[gffGeneIdAttributeName]]),
                        stringsAsFactors=FALSE)
      gns <- gns[gns$type=="mRNA",]
      gns <- gns[,c("tx_name","gene_id")]
      ## Then merge in the more reliable tx_ids
      gns <- merge(gns, transcripts, by="tx_name")
      gns <- gns[,c("tx_id","gene_id")]
    }else{
      warning("No gene information present in gff file")
    }
  }else{
    if(length(gns$ID) != length(unique(gns$ID))){
      stop("Unexpected gene duplicates")}
    ## After testing for genes, I get the actual data from mRNA rows...
    ## The only difference is that in this more normal case the Parents of
    ## these rows will be the genes that I detected previously.
    txs <- .prepareGFF3TXS(data)
    txsGene <- txs[,c("tx_id","Parent")]
    names(txsGene) <- c("tx_id","gene_id")
    ## Then subset by gene_ids
    gns <- txsGene[txsGene$gene_id %in% gns$ID,]
  }
  as.data.frame(gns)
}

.prepareGFF3Fragments <- function(data, type){
  res <- data[data$type==type,]
  if(dim(res)[1] < 1){stop(paste("No",type,"information present in gff file"))
  }else{
    name <- paste0(tolower(type), "_id")
    res <- data.frame(res, data.frame(name=1:dim(res)[1]),
                      stringsAsFactors=FALSE)
    colnames <- c("XXX_chrom","XXX_start","XXX_end","XXX_strand","type",
                    "XXX_name","tx_name","exon_rank","XXX_id")
    names(res) <- sub("XXX", tolower(type), colnames)
  }
  res
}


.prepareGFF3Tables  <- function(gff,exonRankAttributeName,
                                gffGeneIdAttributeName){
  ## pre-clean the data
  data <- .prepareGFF3data.frame(gff,exonRankAttributeName,
                                gffGeneIdAttributeName)
  
  tables <- list()
  ## Get transcripts
  transcripts <- .prepareGFF3transcripts(data)
  tables[[1]] <- transcripts
  names(tables)[1] = "transcripts"
  ## Get genes
  tables[[2]] <- .prepareGFF3genes(data, transcripts,
                                   gffGeneIdAttributeName, gff)
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
  splicings <- splicings[,c('exon_rank','exon_id','exon_name','exon_chrom',
                            'exon_strand','exon_start','exon_end','cds_id',
                            'cds_name','cds_start','cds_end','tx_name')]
  txIds <- unique(transcripts[,c("tx_id","tx_name")])
  splicings <- merge(txIds, splicings, by="tx_name")[,-1]
  ## Clean up any NA columns (any that are optional)
  splicings <- .cleanSplicingsNAs(splicings)  
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
                       type=as.character(gff$type),
                       gene_id=as.character(gff$gene_id),
                       transcript_id=as.character(gff$transcript_id),
                       stringsAsFactors=FALSE)
  ## add ExonRank if there is any  
  if(!is.null(exonRankAttributeName)){
    gff <- as(gff, "GRanges")
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
                                    exonRankAttributeName=NULL,
                                    gffGeneIdAttributeName=NULL,
                                    chrominfo,
                                    dataSource,
                                    species,
                                    circ_seqs=DEFAULT_CIRC_SEQS,
                                    miRBaseBuild=NULL)
{
  ## Some argument checking:
  if(missing(dataSource)) stop("No Datasource provided")
  if(missing(species)) stop("No species provided")
  format <- match.arg(format)
  
  ## start by importing the relevant features from the specified file
  feature.type <- c("gene", "mRNA", "exon", "CDS")
  gff <- import(file, format=format, feature.type=feature.type,
                asRangedData=TRUE)

  if(format=="gff3"){
    ## check that we have ID, Parent
    if(all(c("ID","Parent") %in% colnames(gff))){
      tables <- .prepareGFF3Tables(gff, exonRankAttributeName,
                                   gffGeneIdAttributeName)
      ## results come back in list like: tables$transctripts etc.
    }
  }else if(format=="gtf"){
    ## check that we have gene_id and transcript_id
    if(all(c("gene_id","transcript_id")
           %in% colnames(gff))){
      tables <- .prepareGTFTables(gff,exonRankAttributeName)
    }
  }
  ## TODO: verify that I have all I really need for metadata
  ## build up the metadata
  metadata <- .prepareGFFMetadata(dataSource,species,miRBaseBuild)

  ## If there is not chrominfo, then make one up best you can (no lengths)
  if(missing(chrominfo)){
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
##                                 exonRankAttributeName = "nb_exon",
##                                 circ_seqs=DEFAULT_CIRC_SEQS,
##                                 miRBaseBuild=NULL)
## saveDb(txdb,file="TESTGFF.sqlite")



## ## TESTING GTF
## gtfFile=system.file("extdata","Aedes_aegypti.partial.gtf")
## txdb <- makeTranscriptDbFromGFF(file=gtfFile,
##                                 format="gtf",
##                                 dataSource="ftp://ftp.ensemblgenomes.org/pub/metazoa/release-13/gtf/aedes_aegypti/",
##                                 species="Aedes aegypti",
##                                 circ_seqs=DEFAULT_CIRC_SEQS,
##                                 miRBaseBuild=NULL)
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
##           species="Drosophila melanogaster",
##           circ_seqs=DEFAULT_CIRC_SEQS,
##           miRBaseBuild=NULL)
## foo <- transcriptsBy(txdb3, by="gene")

## this file will not have any proper names etc.  Use it for testing.
## flyFile = "dmel-1000-r5.11.filtered.gff"


