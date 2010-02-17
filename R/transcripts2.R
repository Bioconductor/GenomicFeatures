## convert a named list into a SQL where condition
.sqlWhereClause <- function(vals)
{
    sql <-
      lapply(seq_len(length(vals)), function(i) {
                 v <- vals[[i]]
                 if (!is.numeric(v))
                     v <- paste("'", v, "'", sep="")
                 v <- paste("(", paste(v, collapse=","), ")", sep="")
                 v <- paste(names(vals)[i], " IN ", v, sep="")
                 paste("(", v, ")", sep="")
             })
    clause <- paste(unlist(sql), collapse = " AND ")
    if (nchar(clause) > 0)
        clause <- paste("WHERE", clause)
    clause
}


transcripts2 <- function(txdb, partition=c("none", "exon", "cds"), filter=NULL)
{
    ## check to see if user wanted deprecated function
    if(is.data.frame(txdb))
        stop("Please use 'transcripts_deprecated' for older data.frame-based transcript metadata.")

    ## check that txdb is a TranscriptDb object
    if(!is(txdb,"TranscriptDb"))
        stop("'txdb' must be a TranscriptDb object")

    ## check the filter argument
    validFilterNames <- c("gene_id", "tx_id", "tx_name", "tx_chrom", "tx_strand")
    if(!is.null(filter) &&
       (!is.list(filter) || is.null(names(filter)) ||
        !all(names(filter) %in% validFilterNames))) {
        stop("'filter' must be NULL or a list with names being a combination of ",
             paste(dQuote(validFilterNames), collapse = ", "))
    }
    whichId <- which(names(filter) == "tx_id")
    if(length(whichId) > 0) {
        names(filter)[whichId] <- "t._tx_id"
    }

    partition <- match.arg(partition)
    if (partition == "none") {
        ## create SQL query
        sql <-
          paste("SELECT tx_chrom, tx_start, tx_end, tx_strand,",
                "gene_id, tx_name, t._tx_id AS tx_id",
                "FROM transcript AS t INNER JOIN transcript_rtree AS trt",
                "ON t._tx_id=trt._tx_id",
                "LEFT OUTER JOIN gene AS g",
                "ON t._tx_id=g._tx_id",
                .sqlWhereClause(filter),
                "ORDER BY t._tx_id")

        ## get the data from the database
        ans <- dbGetQuery(txdb@conn, sql)
        ans <-
          split(GenomicFeature(seqnames = ans[["tx_chrom"]],
                               ranges = IRanges(start = ans[["tx_start"]],
                                                end = ans[["tx_end"]]),
                               strand = strand(ans[["tx_strand"]]),
                               gene_id = ans[["gene_id"]],
                               tx_name = ans[["tx_name"]],
                               tx_id = ans[["tx_id"]]),
                ans[["tx_id"]])
    } else if (partition == "exon") {
        ## create SQL query
        sql <-
          paste("SELECT exon_chrom, exon_start, exon_end, exon_strand,",
                "exon_name, e._exon_id AS exon_id, spl._tx_id AS tx_id",
                "FROM splicing AS spl LEFT OUTER JOIN exon AS e",
                "ON spl._exon_id=e._exon_id",
                "LEFT OUTER JOIN exon_rtree AS ert",
                "ON e._exon_id=ert._exon_id",
                "LEFT OUTER JOIN transcript AS t",
                "ON spl._tx_id=t._tx_id",
                "LEFT OUTER JOIN gene AS g",
                "ON spl._tx_id=g._tx_id",
                .sqlWhereClause(filter),
                "ORDER BY spl._tx_id, exon_rank")

        ## get the data from the database
        ans <- dbGetQuery(txdb@conn, sql)
        ans[["tx_id"]] <- factor(ans[["tx_id"]])
        ans <- ans[!is.na(ans[["exon_start"]]),,drop=FALSE]
        ans <-
          split(GenomicFeature(seqnames = ans[["exon_chrom"]],
                               ranges = IRanges(start = ans[["exon_start"]],
                                                end = ans[["exon_end"]]),
                               strand = strand(ans[["exon_strand"]]),
                               exon_name = ans[["exon_name"]],
                               exon_id = ans[["exon_id"]]),
                ans[["tx_id"]])
    } else if (partition == "cds") {
        ## create SQL query
        sql <-
          paste("SELECT cds_chrom, cds_start, cds_end, cds_strand,",
                "cds_name, c._cds_id AS cds_id, spl._tx_id AS tx_id",
                "FROM splicing AS spl LEFT OUTER JOIN cds AS c",
                "ON spl._cds_id=c._cds_id",
                "LEFT OUTER JOIN cds_rtree AS crt",
                "ON c._cds_id=crt._cds_id",
                "LEFT OUTER JOIN transcript AS t",
                "ON spl._tx_id=t._tx_id",
                "LEFT OUTER JOIN gene AS g",
                "ON spl._tx_id=g._tx_id",
                .sqlWhereClause(filter),
                "ORDER BY spl._tx_id, exon_rank")

        ## get the data from the database
        ans <- dbGetQuery(txdb@conn, sql)
        ans[["tx_id"]] <- factor(ans[["tx_id"]])
        ans <- ans[!is.na(ans[["cds_start"]]),,drop=FALSE]
        ans <-
          split(GenomicFeature(seqnames = ans[["cds_chrom"]],
                               ranges = IRanges(start = ans[["cds_start"]],
                                                end = ans[["cds_end"]]),
                               strand = strand(ans[["cds_strand"]]),
                               cds_name = ans[["cds_name"]],
                               cds_id = ans[["cds_id"]]),
                ans[["tx_id"]])
    }
    ans
}
