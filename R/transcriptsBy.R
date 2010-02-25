##  Helper methods
.splitField <- function(groupBy){
  sf <- switch(groupBy,
               "tx"   = "_tx_id",
               "exon" = "_exon_id",
               "cds"  = "_cds_id",
               "gene" = "gene_id")
  sf
}

.checkByArgs <- function(txdb, groupBy){
  ## check that txdb is a TranscriptDb object
  if(!is(txdb,"TranscriptDb"))
    stop("'txdb' must be a TranscriptDb object")

  ## check the groupBy argument
  validGroupBy <- c("tx", "gene", "exon","cds")
  if(length(groupBy) == 1 &&
     (!is.character(groupBy) || !all(groupBy %in% validGroupBy))) {
    stop("'groupBy' must be one of ",
         paste(dQuote(validGroupBy), collapse = ", "))
  }
}


## Lets do the obvious in the most straightforward of ways.
transcriptsBy <- function(txdb, groupBy="gene")
{
  .checkByArgs(txdb, groupBy)
  ## basic transcriptBy query
  sql <- switch(groupBy,
    "tx" = "SELECT * FROM transcript AS t, transcript_rtree AS tr
            WHERE t._tx_id=tr._tx_id ORDER BY t._tx_id",
    "exon" = "SELECT * FROM transcript AS t, transcript_rtree AS tr,
              splicing AS s, exon AS e
              WHERE t._tx_id=tr._tx_id AND t._tx_id=s._tx_id AND
              e._exon_id=s._exon_id ORDER BY e._exon_id",
    "cds" = "SELECT * FROM transcript AS t, transcript_rtree AS tr,
             splicing AS s, cds AS c
             WHERE t._tx_id=tr._tx_id AND t._tx_id=s._tx_id AND
             c._cds_id=s._cds_id ORDER BY c._cds_id",
    "gene" =  "SELECT * FROM transcript AS t, transcript_rtree AS tr,
               gene AS g
               WHERE t._tx_id=tr._tx_id AND t._tx_id=g._tx_id
               ORDER BY g.gene_id"
  )
  
  ## get the data
  ans <- dbGetQuery(txdb@conn, sql)

  ## package the data into a gRangesList,
  ## splitting on the collumn that you grouped by
  splitField <- .splitField(groupBy)
  otherFields <- switch(groupBy,
                        "tx"   = "tx_name",
                        "exon" = c("tx_name", "exon_rank", "exon_name",
                                   "exon_chrom", "exon_strand"),
                        "cds"  = c("tx_name", "exon_rank", "cds_name",
                                   "cds_chrom", "cds_strand"),
                        "gene" = "tx_name")

  ##Now I need to make a GRanges object
  grs <- GRanges(seqnames = ans[["tx_chrom"]],
                 ranges   = IRanges(start = ans[["tx_start"]],
                   end    = ans[["tx_end"]]),
                 strand   = strand(ans[["tx_strand"]]),
                 ans[,otherFields, drop=FALSE])
                 
  ## Then split into a GRangesList (which is what we return)
  split(grs, as.factor(ans[[splitField]]))
}



exonsBy <- function(txdb, groupBy="gene")
{
  .checkByArgs(txdb, groupBy)
  ## basic transcriptBy query
  sql <- switch(groupBy,
    "tx" = "SELECT * FROM exon AS e, exon_rtree AS er,
            splicing AS s, transcript AS t
            WHERE e._exon_id=er._exon_id AND e._exon_id=s._exon_id AND
            t._tx_id=s._tx_id ORDER BY t._tx_id",
    "exon" = "SELECT * FROM exon AS e, exon_rtree AS er
              WHERE e._exon_id=er._exon_id  ORDER BY e._exon_id",
    "cds" = "SELECT * FROM exon AS e, exon_rtree AS er,
             splicing AS s, cds AS c
             WHERE e._exon_id=er._exon_id AND e._exon_id=s._exon_id AND
             c._cds_id=s._cds_id ORDER BY c._cds_id",
    "gene" = "SELECT * FROM exon AS e, exon_rtree AS er,
              splicing AS s, transcript AS t, gene AS g
              WHERE e._exon_id=er._exon_id AND e._exon_id=s._exon_id AND
              t._tx_id=s._tx_id AND g._tx_id=t._tx_id ORDER BY g.gene_id",
  )
  
  ## get the data
  ans <- dbGetQuery(txdb@conn, sql)

  ## package the data into a gRangesList,
  ## splitting on the collumn that you grouped by
  splitField <- .splitField(groupBy)
  otherFields <- switch(groupBy,
                        "tx"   = c("exon_name", "exon_rank", "tx_name",
                                   "tx_chrom", "tx_strand"),
                        "exon" = "exon_name",
                        "cds"  = c("exon_name", "exon_rank", "cds_name",
                                   "cds_chrom", "cds_strand"),
                        "gene" = c("exon_name", "exon_rank", "tx_name",
                                   "tx_chrom", "tx_strand"))

  ##Now I need to make a GRanges object
  grs <- GRanges(seqnames = ans[["exon_chrom"]],
                 ranges   = IRanges(start = ans[["exon_start"]],
                   end    = ans[["exon_end"]]),
                 strand   = strand(ans[["exon_strand"]]),
                 ans[,otherFields, drop=FALSE])
                 
  ## Then split into a GRangesList (which is what we return)
  split(grs, as.factor(ans[[splitField]]))
}





cdsBy <- function(txdb, groupBy="gene")
{
  .checkByArgs(txdb, groupBy)
  ## basic transcriptBy query
  sql <- switch(groupBy,
    "tx" = "SELECT * FROM cds AS c, cds_rtree AS cr,
            splicing AS s, transcript AS t
            WHERE c._cds_id=cr._cds_id AND c._cds_id=s._cds_id AND
            t._tx_id=s._tx_id ORDER BY t._tx_id",
    "exon" = "SELECT * FROM cds AS c, cds_rtree AS cr,
             splicing AS s, exon AS e
             WHERE c._cds_id=cr._cds_id AND c._cds_id=s._cds_id AND
             e._exon_id=s._exon_id ORDER BY e._exon_id",
    "cds" = "SELECT * FROM cds AS c, cds_rtree AS cr
              WHERE c._cds_id=cr._cds_id  ORDER BY c._cds_id",
    "gene" = "SELECT * FROM cds AS c, cds_rtree AS cr,
              splicing AS s, transcript AS t, gene AS g
              WHERE c._cds_id=cr._cds_id AND c._cds_id=s._cds_id AND
              t._tx_id=s._tx_id AND g._tx_id=t._tx_id ORDER BY g.gene_id",
  )
  
  ## get the data
  ans <- dbGetQuery(txdb@conn, sql)

  ## package the data into a gRangesList,
  ## splitting on the collumn that you grouped by
  splitField <- .splitField(groupBy)
  otherFields <- switch(groupBy,
                        "tx"   = c("cds_name", "exon_rank", "tx_name",
                                   "tx_chrom", "tx_strand"),
                        "exon"  = c("cds_name", "exon_rank", "exon_name",
                                   "exon_chrom", "exon_strand"),
                        "cds" = "cds_name",
                        "gene" = c("cds_name", "exon_rank", "tx_name",
                                   "tx_chrom", "tx_strand"))

  ##Now I need to make a GRanges object
  grs <- GRanges(seqnames = ans[["cds_chrom"]],
                 ranges   = IRanges(start = ans[["cds_start"]],
                   end    = ans[["cds_end"]]),
                 strand   = strand(ans[["cds_strand"]]),
                 ans[,otherFields, drop=FALSE])
                 
  ## Then split into a GRangesList (which is what we return)
  split(grs, as.factor(ans[[splitField]]))
}




##tests
## library(GenomicFeatures)
## txdb = loadFeatures("UCSC_knownGene_sample.sqlite")
## groupBy = "exon"
## foo = transcriptsBy(txdb, groupBy)
##works for all but "by gene"
## bar = exonsBy(txdb, groupBy) 
## sna = cdsBy(txdb, groupBy)

## TODO: add unit tests for these, but ALSO do some sanity checking.
