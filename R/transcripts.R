
##Time for methods to access data via more direct queries.

##1st some helper methods for complex restrictions in where clauses.
.orAndTester <- function(val){  
  if(length(unlist(val))>1){
    .appendValElementsSQL(val)
  }else{
    paste("AND ", names(val), " = ","'",val,"'",sep="")
  }
}

.appendValsSQL <- function(vals){
  sql <- character()
  for(i in seq_len(length(vals))){
    #print(i)
    sql <- c(sql, .orAndTester(vals[i]))    
  }
  paste(sql, collapse=" ")
}

.appendValElementsSQL <- function(vals){
  vals = unlist2(vals) ##AnnotationDbi:::unlist2 preserves names with repeats
  sql <- character()
  for(i in seq_len(length(vals))){
    sql <- c(sql, paste(names(vals[i])," = ","'",vals[i],"'",sep=""))
  }
  sql <- paste(sql, collapse=" OR ")
  paste("AND (",sql,")",collapse=" ")
}


## This is the core function for looking up transcripts

transcripts <- function(txdb, vals, columns=c("tx_id", "tx_name"))
{
  ##Add error here
  if(is.data.frame(txdb) && is.integer(vals)){stop("This is not the transcripts function that you are looking for. Please use transcripts_deprecated instead.")}

  ## check that txdb is in fact a TranscriptDb object
  if(!is(txdb,"TranscriptDb"))stop("txdb MUST be a TranscriptDb object.")
  
  ## check the vals:
  valNames <- c("gene_id", "tx_id", "tx_name", "tx_chrom", "tx_strand")
  if(!all(names(vals) %in% valNames)){
    stop(paste("Argument names for vals must be some combination of: ",
               valNames,sep=""))
  }

  ## check the cols:
  colNames <- c("tx_id", "tx_name", "gene_id", "exon_id","cds_id",
                NULL, character(0))
  if(!all(columns %in% colNames)){
    stop(paste("Arguments to column must be some combination of: ",
               colNames,sep=""))
  }

  ## base SQL query 1
  sql <- paste("SELECT gene_id, t._tx_id AS tx_id, tx_name, tx_chrom, tx_strand,",
               "tx_start, tx_end, _exon_id AS exon_id, _cds_id AS cds_id",
               "FROM transcript AS t, transcript_rtree AS trt,",
               "gene AS g, splicing AS s",
               "WHERE t._tx_id=trt._tx_id AND t._tx_id=g._tx_id",
               "AND t._tx_id=s._tx_id")


  ## Now we just need to finish the where clause with stuff in "vals"
  sql <- paste(sql, .appendValsSQL(vals), collapse=" ")
  ans <- dbGetQuery(txdb@conn, sql)
  
  if(dim(ans)[1] >0){
      rd <- .formatRD(ans, "get", "tx")
      if(is.null(columns) || length(columns)==0){
        return(rd)
      }else{
        return(.appendCols(rd, ans, columns))
      }
  }else{warning("Please be advised that no matching data was found.")}
}






## This is the core function for looking up exons

exons <- function(txdb, vals)
{
  ##Add error here
  if(is.data.frame(txdb) && is.integer(vals)){stop("This is not the exons function that you are looking for. Please use exons_deprecated instead.")}

  ## check that txdb is in fact a TranscriptDb object
  if(!is(txdb,"TranscriptDb"))stop("txdb MUST be a TranscriptDb object.")
  
  ## check the vals:
  valNames <- c("exon_id", "exon_strand", "tx_chrom", "tx_strand",
                "tx_id", "tx_name")
  if(!all(names(vals) %in% valNames)){
    stop(paste("Argument names for vals must be some combination of: ",
               valNames,sep=""))
  }

  ## base SQL query 1  -- JOIN THIS to transcript
  sql <- paste("SELECT e._exon_id AS exon_id, exon_chrom, exon_strand,",
               "exon_start, exon_end, tx_chrom, tx_strand, tx_name,",
               "t._tx_id AS tx_id",
               "FROM exon AS e, exon_rtree AS ert, transcript AS t,",
               "splicing AS s",
               "WHERE e._exon_id=ert._exon_id AND e._exon_id=s._exon_id",
               "AND t._tx_id=s._tx_id")

  ## Now we just need to finish the where clause with stuff in "vals"
  sql <- paste(sql, .appendValsSQL(vals), collapse=" ")
  ans <- dbGetQuery(txdb@conn, sql)

  ## We always return the exon_id
  columns <- "exon_id"
  if(dim(ans)[1] >0){
      rd <- .formatRD(ans, "get", "exon")
      return(.appendCols(rd, ans, columns))
  }else{warning("Please be advised that no matching data was found.")}
}


