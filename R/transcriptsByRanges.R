## Get the list of possible values for a type from a table:
.getValsFromTable <- function(txann, colname, table){
  sql <- paste("SELECT", colname, "FROM", table)
  as.character(unlist(unique(dbGetQuery(txann@conn, sql))))
}

## Check that the type of thing (chromosome, strand etc.) is in the
## transcript table
.checkFields <- function(txann, colname, data, table){
  annot <- .getValsFromTable(txann, colname, table) 
  all(data %in% annot)
}


## in the absence of chromosome and strand information, expand the query...
.createUniqueChromStrand <- function(txdb, table, chrom, strand, colnames){
  if(is.null(chrom) && !is.null(strand)){
    allChrms = .getValsFromTable(txdb, colnames[1L], table)
    uChromStrand <- unique(cbind(rep(allChrms,length(strand)),
                                 rep(strand,each=length(allChrms))) )  
  }
  if(!is.null(chrom) && is.null(strand)){
    allStrnds = .getValsFromTable(txdb, colnames[2L], table)
    uChromStrand <- unique(cbind(rep(chrom,length(allStrnds)),
                                 rep(allStrnds,each=length(chrom))) )  
  }
  if(is.null(strand) && is.null(chrom)){
    allChrms = .getValsFromTable(txdb, colnames[1L], table)
    allStrnds = .getValsFromTable(txdb, colnames[2L], table)
    uChromStrand <- unique(cbind(rep(allChrms,length(allStrnds)),
                                 rep(allStrnds,each=length(allChrms))))
  }
  if(!is.null(chrom) && !is.null(strand)){
    uChromStrand <- unique(cbind(as.character(chrom),as.character(strand)))
  }
  uChromStrand
}


.setRange <- function(ranges, txdb, chrom, strand, chr, strnd){
  if(is.null(chrom) && !is.null(strand)){
    range <- ranges[strand==strnd]
  }
  if(!is.null(chrom) && is.null(strand)){
    range <- ranges[chrom==chr]
  }
  if(is.null(strand) && is.null(chrom)){
    range <- ranges
  }
  if(!is.null(chrom) && !is.null(strand)){
    range <- ranges[chrom==chr & strand==strnd]
  }
  range
}


.mapFeatures <- function(txdb, sql, range, chrom, strnd, colprefix,
                         type="any"){

  sqlChrom <- paste(" AND ",colprefix,"_chrom = '",chrom,"'",sep="")
  sqlStrand <- paste(" AND ",colprefix,"_strand = '",strnd,"'",sep="")
  sql2 <- paste(sql, sqlChrom, sqlStrand)
  res <- dbGetQuery(txdb@conn, sql2)
  
  ## Ranges is the query, the DB returned tx_start and tx_end will be the
  ## subject
  startName <- paste(colprefix,"_start",sep="")
  endName <- paste(colprefix,"_end",sep="")
  if(dim(res)[1] > 0){
    map <- findOverlaps(range,
                        IRanges(start = unlist(res[startName]),
                                end = unlist(res[endName])),
                        multiple=TRUE,
                        type=type)
  
    ##Then subset the query result with only the matching stuff.
    ##And I also have to match that up with the stuff in range...
    matchRange <- range[map@matchMatrix[,1]]
    if(!is.na(map@matchMatrix[1])){
      ans <- cbind(as.data.frame(matchRange),res[map@matchMatrix[,2],])
    }else{ans <- data.frame()}
  }else{ans <- data.frame()}
  ans
}


##Format the ranged data object depending on whether mapping or getting
.formatRD <- function(ans, format, prefix){
  switch(format,
         "map" = return(RangedData(
           ranges = IRanges(start = ans[["start"]],
                              end = ans[["end"]]),
           strand = factor(ans[[paste(prefix,"_strand",sep="")]],
                           levels=c("-","+","*")),
           space = ans[[paste(prefix,"_chrom",sep="")]])),
         "get" = return(RangedData(
           ranges = IRanges(start = ans[[paste(prefix,"_start",sep="")]],
                              end = ans[[paste(prefix,"_end",sep="")]]),
           strand = factor(ans[[paste(prefix,"_strand",sep="")]],
                           levels=c("-","+","*")),
           space = ans[[paste(prefix,"_chrom",sep="")]]))         
         )
}

##Add addional data that was requested to the Ranged Data that is returned.
.appendCols <- function(rd, ans, col){## for each value in cols
  ## verify that cols are in ans
  if(all(col %in% colnames(ans))){
    for(i in seq_len(length(col))){
      rd[[col[i]]] <- ans[[col[i]]]
    }
  }else{stop(paste("There is either a problem with the requested columns",
                   "or the answer returned from the DB is not complete."))}
  rd
}


transcriptsByRanges <- function(txdb, ranges, restrict = "any",
                                columns=c("tx_id", "tx_name"))
{
  ## check/assign the strand, ranges c-somes
  if(.checkFields(txdb, "tx_strand", ranges[["strand"]], "transcript")){
    strand <- ranges[["strand"]]
  }else{stop("Strand values for ranges do not match annotation DB")}
  if(.checkFields(txdb, "tx_chrom", space(ranges), "transcript")){
    chrom <- space(ranges)
  }else{stop("Space values for ranges do not match annotation DB")}
  ranges <- unlist(ranges(ranges), use.names=FALSE)
  
  ## check that txdb is in fact a TranscriptDb object
  if(!is(txdb,"TranscriptDb"))stop("txdb MUST be a TranscriptDb object.")
  ## check the columnss:
  colNames <- c("tx_id", "tx_name", "gene_id", "exon_id","cds_id",
                NULL, character(0))
  if(!all(columns %in% colNames)){
    stop(paste("Arguments to column must be some combination of: ",
               colNames,sep=""))
  }
  ## base SQL query
  sql <- paste("SELECT gene_id, t._tx_id AS tx_id, tx_name, tx_chrom, tx_strand,",
               "tx_start, tx_end, _exon_id AS exon_id, _cds_id AS cds_id",
               "FROM transcript AS t, transcript_rtree AS trt,",
               "gene AS g, splicing AS s",
               "WHERE t._tx_id=trt._tx_id AND t._tx_id=g._tx_id",
               "AND t._tx_id=s._tx_id")
  ## We could speed this up by doing less complicated joins when possible...
  
  uChromStrand <- .createUniqueChromStrand(txdb,"transcript",
                                           chrom, strand,
                                           c("tx_chrom", "tx_strand"))
  ans <- data.frame()
  for(i in seq_len(dim(uChromStrand)[1])){
    chr <- uChromStrand[i,1]
    strnd <- uChromStrand[i,2]
    range <- .setRange(ranges, txdb, chrom, strand, chr, strnd)
    if(length(range) > 0){
      newData <- .mapFeatures(txdb, sql, range, chr, strnd, "tx", restrict)
      if(dim(newData)[1] > 0){
        chrs <- rep(chr, dim(newData)[1])
        strnds <-rep(strnd, dim(newData)[1])
        ans <- rbind(ans, cbind(chrs,strnds,newData))
      }
    }
  }
  if(dim(ans)[1] >0){
      rd <- .formatRD(ans, "get", "tx")
      if(is.null(columns) || length(columns)==0){
        return(rd)
      }else{
        return(.appendCols(rd, ans, columns))
      }
  }else{warning("Please be advised that no matching data was found.")}
}











##Additional code to support mapping on exons:

## TODO: Need to handle the case where we don't have a ranges... (missing
## method is gone)
exonsByRanges <- function(txdb, ranges, restrict = "any",
                          columns=c("exon_id", "exon_name"))
{
  ## check/assign the strand, ranges c-somes
  if(.checkFields(txdb, "exon_strand", ranges[["exon_strand"]], "exon")){
    strand <- ranges[["strand"]]
  }else{stop("Strand values for ranges do not match annotation DB")}
  if(.checkFields(txdb, "exon_chrom", space(ranges), "exon")){
    chrom <- space(ranges)
  }else{stop("Space values for ranges do not match annotation DB")}
  ranges <- unlist(ranges(ranges), use.names=FALSE)
  
  ## check that txdb is in fact a TranscriptDb object
  if(!is(txdb,"TranscriptDb"))stop("txdb MUST be a TranscriptDb object.")

  ## check the columnss:
  colNames=c("exon_id", NULL, character(0))
  if(!all(columns %in% colNames)){
    stop(paste("Arguments to column must be some combination of: ",
               colNames,sep=""))
  }  
  ## base SQL query:
  sql <- paste("SELECT exon_id, exon_chrom, exon_strand,",
               "exon_start, exon_end",
               "FROM exon INNER JOIN exon_rtree",
               "ON (exon._exon_id=exon_rtree._exon_id")
  
  uChromStrand <- .createUniqueChromStrand(txdb,"exon",
                                           chrom, strand,
                                           c("exon_chrom", "exon_strand"))
  ans <- data.frame()
  for(i in seq_len(dim(uChromStrand)[1])){
    chr <- uChromStrand[i,1]
    strnd <- uChromStrand[i,2]
    range <- .setRange(ranges, txdb, chrom, strand, chr, strnd)
    if(length(range) > 0){
      newData <- .mapFeatures(txdb, sql, range, chr, strnd , "exon", restrict)
      if(dim(newData)[1] > 0){
        chrs <- rep(chr, dim(newData)[1])
        strnds <-rep(strnd, dim(newData)[1])
        ans <- rbind(ans, cbind(chrs,strnds,newData))
      }
    }
  }
  if(dim(ans)[1] >0){
      rd <- .formatRD(ans, "get", "exon")
      if(is.null(columns) || length(columns)==0){
        return(rd)
      }else{
        return(.appendCols(rd, ans, columns))
      }    
    }else{warning("Please be advised that no matching data was found.")}
}



