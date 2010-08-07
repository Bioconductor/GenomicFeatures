## convert a named list into an SQL where condition
.sqlWhereIn <- function(vals)
{
  if (length(vals) == 0) {
    character(0)
  } else {
    sql <-
      lapply(seq_len(length(vals)), function(i) {
               v <- vals[[i]]
               if (!is.numeric(v))
                 v <- paste("'", v, "'", sep="")
               v <- paste("(", paste(v, collapse=","), ")", sep="")
               v <- paste(names(vals)[i], " IN ", v, sep="")
               paste("(", v, ")", sep="")
            })
    paste("WHERE", paste(unlist(sql), collapse = " AND "))
  }
}


## transcripts function and helpers

.newListBySplit <- function(class, x, f)
{
  IRanges:::newCompressedList(class, unlistData = x, splitFactor = f)
}

.geneCharacterList <- function(txdb, tx_ids)
{
  sqlIDs <- paste("(", paste(tx_ids, collapse=","), ")", sep="")
  sql <- paste("SELECT gene_id, _tx_id FROM gene WHERE _tx_id IN ", sqlIDs,
               sep="")
  ans <- dbEasyQuery(txdbConn(txdb), sql)
  ans[["_tx_id"]] <-
    factor(as.character(ans[["_tx_id"]]), levels=as.character(tx_ids))
  ans <- ans[order(ans[["_tx_id"]]), ,drop=FALSE]
  unname(.newListBySplit("CompressedCharacterList", ans[["gene_id"]], ans[["_tx_id"]]))
}

.exonORcdsIntegerList <- function(txdb, tx_ids, type=c("exon", "cds"))
{
  type <- match.arg(type)
  type_id <- paste("_", type, "_id", sep="")
  sqlIDs <- paste("(", paste(tx_ids, collapse=","), ")", sep="")
  sql <- paste("SELECT _tx_id, exon_rank, ", type_id, " ",
               "FROM splicing WHERE _tx_id IN ", sqlIDs, sep="")
  ans <- dbEasyQuery(txdbConn(txdb), sql)
  ans[["_tx_id"]] <-
    factor(as.character(ans[["_tx_id"]]), levels=as.character(tx_ids))
  ans <- ans[order(ans[["_tx_id"]], ans[["exon_rank"]]), ,drop=FALSE]
  .newListBySplit("CompressedIntegerList", ans[[type_id]], ans[["_tx_id"]])
}


transcripts <- function(txdb, vals=NULL, columns=c("tx_id", "tx_name"))
{
  ## check to see if user wanted deprecated function
  if(is.data.frame(txdb))
    stop("Please use 'transcripts_deprecated' for older data.frame-based transcript metadata.")

  ## check that txdb is a TranscriptDb object
  if(!is(txdb,"TranscriptDb"))
    stop("'txdb' must be a TranscriptDb object")

  ## check the vals argument
  validValNames <- c("gene_id", "tx_id", "tx_name", "tx_chrom", "tx_strand")
  if(!is.null(vals) &&
     (!is.list(vals) || is.null(names(vals)) ||
      !all(names(vals) %in% validValNames))) {
    stop("'vals' must be NULL or a list with names being a combination of ",
         paste(dQuote(validValNames), collapse = ", "))
  }
  whichId <- which(names(vals) == "tx_id")
  if(length(whichId) > 0) {
    names(vals)[whichId] <- "transcript._tx_id"
  }

  ## check the columns argument
  validColumns <- c("tx_id", "tx_name", "gene_id", "exon_id","cds_id")
  if(length(columns) > 0 &&
     (!is.character(columns) || !all(columns %in% validColumns))) {
    stop("'columns' must be NULL or a combination of ",
         paste(dQuote(validColumns), collapse = ", "))
  }

  ## create SQL query
  if ("tx_name" %in% columns)
    optionalColumn <- ", tx_name"
  else
    optionalColumn <- ""
  if ("gene_id" %in% names(vals))
    optionalLeftJoin <- "LEFT JOIN gene ON transcript._tx_id=gene._tx_id"
  else
    optionalLeftJoin <- ""
  sql <- paste("SELECT tx_chrom, tx_start, tx_end, tx_strand,",
               "transcript._tx_id AS tx_id", optionalColumn,
               "FROM transcript", optionalLeftJoin,
               .sqlWhereIn(vals),
               "ORDER BY tx_chrom, tx_strand, tx_start, tx_end")

  ## get the data from the database
  ans <- dbEasyQuery(txdbConn(txdb), sql)
  seqlengths <- seqlengths(txdb)
  ans <-
    GRanges(seqnames = factor(ans[["tx_chrom"]], levels = names(seqlengths)),
            ranges = IRanges(start = ans[["tx_start"]],
                             end = ans[["tx_end"]]),
            strand = strand(ans[["tx_strand"]]),
            ans[-c(1:4)],
            seqlengths = seqlengths)

  if(length(ans) > 0 && any(c("gene_id", "exon_id","cds_id") %in% columns)) {
    if("gene_id" %in% columns) {
      elementMetadata(ans)[["gene_id"]] <-
        .geneCharacterList(txdb, elementMetadata(ans)[["tx_id"]])
    }
    if("exon_id" %in% columns) {
      elementMetadata(ans)[["exon_id"]] <-
        .exonORcdsIntegerList(txdb, elementMetadata(ans)[["tx_id"]], "exon")
    }
    if("cds_id" %in% columns) {
      elementMetadata(ans)[["cds_id"]] <-
        .exonORcdsIntegerList(txdb, elementMetadata(ans)[["tx_id"]], "cds")
    }
  }

  if (!("tx_id" %in% columns))
    elementMetadata(ans)[["tx_id"]] <- NULL

  ans
}


## exon and cds functions and helper

.exonORcdsGRanges <- function(txdb, vals=NULL, type=c("exon", "cds"))
{
  type <- match.arg(type)

  if (type == "exon" && is.data.frame(txdb))
    stop("Please use 'exons_deprecated' for older data.frame-based transcript metadata.")

  ## check that txdb is a TranscriptDb object
  if(!is(txdb,"TranscriptDb"))
    stop("'txdb' must be a TranscriptDb object")

  ## check the vals argument
  validValNames <- gsub("TYPE", type, c("TYPE_id", "TYPE_chrom", "TYPE_strand"))
  if(!is.null(vals) &&
     (!is.list(vals) || is.null(names(vals)) ||
      !all(names(vals) %in% validValNames))) {
    stop("'vals' must be NULL or a list with names being a combination of ",
         paste(dQuote(validValNames), collapse = ", "))
  }
  whichId <- which(names(vals) == paste(type, "_id", sep=""))
  if(length(whichId) > 0) {
      names(vals)[whichId] <- paste("TYPE._", type, "_id", sep="")
  }

  ## create base SQL query
  sql <- paste("SELECT TYPE_chrom, TYPE_start, TYPE_end, TYPE_strand,",
               "TYPE._TYPE_id AS TYPE_id FROM TYPE",
               .sqlWhereIn(vals),
               "ORDER BY TYPE_chrom, TYPE_strand, TYPE_start, TYPE_end")
  sql <- gsub("TYPE", type, sql)

  ## get the data from the database
  ans <- dbEasyQuery(txdbConn(txdb), sql)
  seqlengths <- seqlengths(txdb)
  ans <-
    GRanges(seqnames =
            factor(ans[[paste(type, "_chrom", sep="")]],
                   levels = names(seqlengths)),
            ranges = IRanges(start = ans[[paste(type, "_start", sep="")]],
                             end = ans[[paste(type, "_end", sep="")]]),
            strand = strand(ans[[paste(type, "_strand", sep="")]]),
            "TYPE_id" = ans[[paste(type, "_id", sep="")]],
            seqlengths = seqlengths)
  colnames(elementMetadata(ans)) <-
    gsub("TYPE", type, colnames(elementMetadata(ans)))

  ans
}


exons <- function(txdb, vals=NULL)
{
  .exonORcdsGRanges(txdb, vals=vals, type="exon")
}


cds <- function(txdb, vals=NULL)
{
  .exonORcdsGRanges(txdb, vals=vals, type="cds")
}
