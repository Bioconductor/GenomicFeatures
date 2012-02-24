## Provide tools for generating a splicing graph

uunlist <- function(x) unlist(x, use.names=FALSE)
ulapply <- function(...) uunlist(lapply(...))
as.char <- function(x) as.character(uunlist(x))

getGH <- function(txdb){
    con <- AnnotationDbi:::dbConn(txdb)
    dbGetQuery(con, "SELECT gene_id from gene")[, 1, drop=TRUE]
}

fastDisjoinOnGRangesList <- function(grl) {
  tmp_gr <- GenomicRanges:::deconstructGRLintoGR(grl)
  tmp_gr2 <- disjoin(tmp_gr)
  GenomicRanges:::reconstructGRLfromGR(tmp_gr2, grl)
}
                                                                      
txclean <- function(tx, seqlevels=paste("chr", 1:21, sep=""), ...){
    ## remove unwanted seqlevels, transcrips spanning multiple
    ## chromosomes / strand
    
    lenunique1 <- function(lst) {
        elts <- uunlist(lst)
        
        ## as.integer(elts) gives pos of objects in the list
        ## if one transcript has the same pos it is removed
        d <- c(0L, cumsum(diff(as.integer(elts)) != 0))
        b <- cumsum(elementLengths(lst))
        d[c(1L, b[-length(b)] + 1L)] == d[b]
    }
    
    ## standard chromosomes
    ## tx <- keepSeqlevels(tx, seqlevels)
    
    ## span single sequence & strand
    seq_ok <- lenunique1(seqnames(tx))

    ## sorts out genes mapping to two different strands
    str_ok <- lenunique1(strand(tx))
    tx <- tx[seq_ok & str_ok]

    ## store info how many genes got removed because of mapping to
    ## different genomic locations
    metadata(tx)[["cleaned"]] <-
      c(seqlevels=sum(!seq_ok), strand=sum(!str_ok))
    tx
}

iodeg <- function(v, e) {
    ## # Add in-degree and out-degree
    m <- max(e)
    
    ## assign in degree to the vertices
    v[["Indeg"]] <- tabulate(e$To, m)[v$Vid]
    
    ## assign out degree to the vertices
    v[["Outdeg"]] <- tabulate(e$From, m)[v$Vid]
    v
}


vertexTranscripts <- function(v) {

    ## get a bool index of vertices which are not root or leaf
    idx <- !v$Type %in% c("0", "1")
    
    ## '0', '1'types get all Tx in Gn
    f <- factor(v$Gn[idx], levels=unique(v$Gn[idx]))
        
    ## list of unique transcript ids per gene
    tscript <- as.vector(sapply(split(v$Tx[idx], f), unique))
    
    ## Give the list of transcripts per genes as name the root and leaf
    ## vertices idex. Assign to the individual vertices the
    ## corresponding transcript ids
    test <- c(setNames(tscript, v$Vid[v$Type=="0"]),
              setNames(tscript, v$Vid[v$Type=="1"]),
              split(v$Tx[idx], v$Vid[idx]))
    
    ## modification
    test[!is.na(names(test))]
}

edgeTranscripts <- function(vtscript, e, e2) {

    ue <- e[e$From %in% e2$From,]
  
    ## extract for each remaining edge those transcripts which
    ## tx id overlap of "from" node and "to" node of the edges
    etscript <- Map(ue$intersect, vtscript[as.character(ue$From)],
                    vtscript[as.character(ue$To)])
        
    ## name the list of related transcripts by the edge names 
    setNames(etscript, ue$Eid)

}


collapseEdge <- function(v, e){
    ## collapse edge with Indeg == Outdeg == 1 
    to <- e$To
    
    ## All vertices with outdegree or indegree equal to 1,
    ## are not interesting since the provide no information about
    ## alternative splicing
    
    ## Root and leave are retained
    ## end ... vertex ids with either out degreen or indegree not equal
    ## to one
    
    end <- unique(v$Vid[v$Indeg!=1 | v$Outdeg != 1])
    
    
    ## start ... vertices with out degree equal 1
    start <- unique(v$Vid[v$Outdeg==1])
    
    ## start not in end which means indegree larger than 1
    ## x contains vertices which are uninformative
    x <- start[!start %in% end]
    
    ## while there are vertices in the graphs with outdegree and indegree
    ## == 1
    
    while (length(x)) {
        i <- match(x, to)
        j <- match(x, e$From)
        
        to[i] <- x <- to[j]
        x <- x[!x %in% end]
    }
    
    i <- e$From %in% c(v$Vid[v$Outdeg!=1], v$Vid[v$Indeg!=1])
    
    e$To <- to
    e <- e[i,]
}


collapseVertex <- function(v) {
    v[!duplicated(v$Vid) & !(v$Indeg == 1 & v$Outdeg == 1),]
}

disjoinExonsByGene <- function(extx, txgn) {
    ## extx exon per transcript
    ## txgn transcripts per gene
    ## group exons into genes, making disjoint
    
    txid <- values(unlist(txgn, use.names=FALSE))[["tx_id"]]
    gnid <- rep(names(txgn), elementLengths(txgn))

    ## names are the transcript ids, content redundant gene ids
    map <- setNames(gnid, txid)

    ## map contains redundant gene IDs and is named by transcript ids
    ## next genes to exons
    gn <- map[rep(names(extx), elementLengths(extx))]

    ## create a exons per gene GRanges list
    exgn <- split(unlist(extx), factor(gn))

    ## Herve gave me that hint for doing the dissjoining much faster
    exgn.orig <- uunlist(exgn)

    ## slow verison
    ## disjoin the exons per gene  
    ## exgn <- endoapply(exgn, disjoin)

    exgn <- fastDisjoinOnGRangesList(exgn)
    
    exgn0 <- unlist(exgn, use.names=FALSE)
    values(exgn0)[["disJ_exon_id"]] <- seq_along(exgn0)

    ## keep the original exon ids after disjoining
    ft <- findOverlaps( exgn.orig, exgn0, select="all", maxgap = 0L)
    orig <-
      split(as.character(values(exgn.orig)[["exon_id"]][queryHits(ft)]),
            values(exgn0)[["disJ_exon_id"]][subjectHits(ft)])
    orig <- orig[values(exgn0)[["disJ_exon_id"]]]
   
    values(exgn0)[["exon_ids"]] <- CharacterList(orig)
    relist(exgn0, exgn)
  }

disjoinExtx <- function(extx, exgn, gTM) {
    ## map exons of extx to (disjoint) exons of exgn
    
    exgn0 <- unlist(exgn, use.names=FALSE)
    extx0 <- unlist(extx, use.names=FALSE)
    values(extx0)[["tx_id"]] <- rep(names(extx), elementLengths(extx))

    gid.t <- gTM[rep(names(extx), elementLengths(extx))]
    gid.g <- rep(names(exgn), elementLengths(exgn))

    olap <- findOverlaps(extx0, exgn0)
    mm <- as.matrix(olap)
    bool <- gid.t[mm[,1]] == gid.g[mm[,2]]
    olap <- olap[bool]
    

    txid <- values(extx0)[["tx_id"]][queryHits(olap)]
    o <- order(subjectHits(olap))
    extx <- split(exgn0[subjectHits(olap)][o], txid[o])

    extx0 <- unlist(extx, use.names=FALSE)
    strand <- unlist(unique(strand(extx)), use.names=FALSE)
    values(extx0)[["exon_rank"]] <- unlist(Map(function(n, s) {
        if (s=="-") rev(seq_len(n)) else seq_len(n)
    }, elementLengths(extx), strand), use.names=FALSE)

    relist(extx0, extx)
}

collapseEdgeExons <- function(v, e) {
    ## group exons by collapsing edges with Indeg == Outdeg == 1
    
    edge <- e$Eid
    to <- e$To

    end <- unique(v$Vid[v$Indeg!=1 | v$Outdeg != 1])
    
    start <- unique(v$Vid[v$Outdeg==1])

    ## same as collpase edges
    x <- start[!start %in% end]

    test <- cbind(edge)
    testT <- cbind(to)
    
    while (length(x)) {
        i <- match(x, to)
        j <- match(x, e$From)
        to[i] <- x <- to[j]
        
        edge[j] <- edge[i]
        x <- x[!x %in% end]
        test <- cbind(test, edge)
        testT <- cbind(testT, to)
    }

    ## to conatins the new ends of the exons
    idx <- match(e$From, v$Vid)
    
    ok <- v$Type[idx] %in% c("[", "-")
    split(v$Ex[idx][ok], edge[ok])
    
}

exonsByEdge <- function(exgn, eexons) {
    ## retrieves the exons assiociated with the individual edges
    ex <- unlist(exgn, use.names=FALSE)
    idx <- match(unlist(eexons, use.names=FALSE),
                 values(ex)$disJ_exon_id)
    
    exsByEdges <-
      split(ex[idx], rep(names(eexons), sapply(eexons, length)))
    ## subset the seqlevels
    seqlevels(exsByEdges) <- as.character(unique(seqnames(ex)))
    exsByEdges
}
    
## creates a data frame of vertices
vertex <- function(extx, txgn, map) {
    
    ttx <- values(uunlist(txgn))[["tx_id"]]
    txnm <- rep(names(txgn), elementLengths(txgn))
    tgn <- setNames(txnm, ttx)
    
    ## Get all genomic start coordinates for the individual genes
    gpos <- local({
        ## index the last transcripts of the individual genes in the 
        ## flat vector
        idx <- cumsum(elementLengths(txgn))
        
        ## get the strand for the transcript with the last index for 
        ## each gene
        strnd <- as.vector(uunlist(strand(txgn))[idx])
        idx <- ifelse(strnd == "+", 1L, -1L)
        mn <- min(start(txgn))
        mx <- max(end(txgn))
        
        ## check if transcript is on the forward or on the reverse
        ## strand
        tst <- idx == 1L
        
        ## get for each gene the minimal/maximal transcript position ->
        ## start (1/-1)
        ## get for each gene the maximal/minimal transcript position ->
        ## end (1/-1)
        
        c(ifelse(tst, mn, mx), ifelse(tst, mx, mn)) * idx
    })

    ## vector of two times the length of the number of genes serving as
    ## root and leave vertices vector: roots <- 0, leaves <- N
    gtx <- c(rep("0", length(txgn)), rep("N", length(txgn)))

    ## initialize type of splice site
    gtype <- factor(c(rep("0", length(txgn)), rep("1", length(txgn))),
                    levels=names(map))

    
    ## Crete a gene data frame providing the most 5' pos of all
    ## transcripts of a certain gene, as well as the gene name, the the
    ## initialized number of transcripts, the initialized number of
    ## exons and the initialized types of the splice sites.
    
    gndf <- data.frame(Pos=gpos, Gn=names(txgn), Tx=gtx, Ex=gtx,
                       Type=gtype, stringsAsFactors=FALSE,
                       keep=rep("K", length(gpos) ))
    
    ## Convert exons per transcript GRanges list object into a flat
    ## GRanges object
    ex <- uunlist(extx)
    

    ## get exon start and exon end positon on both strands
    ## same procedure as with the transcripts
    epos <- local({
        idx <- as.vector(ifelse(strand(ex) == "+", 1L, -1L))
        tst <- idx == 1L
        c(ifelse(tst, start(ex), end(ex)),
            ifelse(tst, end(ex), start(ex))) * idx
    })

    ## get the corresponding exon ids
    eex <- as.character(values(ex)[["disJ_exon_id"]])

    ## repeate the transcript id for all exons of a certain transcript
    ## to create a flat vector of transcript ids corresponding to the
    ## individual exons
    
    etx <- rep(names(extx), elementLengths(extx))


    ## Evaluate the exon type
    ## order of the exons according to their genomic rank
    ## assign the type of the site 
    etype <- local({
        
        ## list of vectors containing the exon ranks per transcript
        rnk <- split(values(ex)[["exon_rank"]],
                     factor(etx, levels=names(extx)))
        
        ## get the bool index of the most 5' exon for each transcript
        first <- values(ex)[["exon_rank"]] == 1
        
        ## get indices of the highest exon ranks for each transcript
        last <- uunlist(lapply(rnk, function(elt) elt == max(elt)))
        
        ## asign the site types to the exon borders      
        factor(c(ifelse(first, "[", "-"), ifelse(last, "]", "^")),
               levels=names(map))
        
    })

    ex.names <- values(ex)[["disJ_exon_id"]]
    uni.idx <- ! duplicated(ex.names)    
    ex.uni <- ex[uni.idx] 
    ex.names <- ex.names[uni.idx]
    
    ex.widths <- width(ranges(ex.uni))
    names(ex.widths) <- ex.names
    
    ## create a exon data frame
    exdf <- data.frame(Pos=epos, Gn=tgn[etx], Tx=etx,
                       Ex=eex, Type=etype, stringsAsFactors=FALSE)
    

    ## check for one NT long exons, if one vertices is assiciated with
    ## such an exon keep the vertices anyway
    leng <- ex.widths[exdf$Ex]
    z <- 1:nrow(exdf)
    k <- rep("K", nrow(exdf))
    
    exdf$keep <- ifelse(unname(ex.widths[exdf$Ex]) > 1, k, z)
    
    
    ## all vertex, ordered by gene, transcript, position, type
    ## combine exon and gene data frame to one single data frame
    ## this means the previous prepared root and leaves get added to
    ## the data frame
    
    df <- rbind(exdf, gndf)
   
    ## sort the data frame first according to the gene id than
    ## according to the transcript id then according to the genomic
    ## position, and then according to the site type factor
    ## 0 (priority 1),[ (priority 2), ^ (priority 3), - (priority 4), ]
    ## (priority 5), 1 (priority 6)

    df <- df[order(df$Gn, df$Pos, df$Type),]
    
    id <- paste(df$Pos, df$Gn, df$Ex, df$keep, sep=":")

    ## check for unique ids --> create a boolean vector
    u <- !duplicated(id)

    ## create a mapping for genearting unique vertex ids by
    ## discarding redundant information since we want an unique
    ## vertex id
    idmap <- setNames(seq_len(sum(u)), id[u])

    ## save vertex ids in the data frame
    df[["Vid"]] <- idmap[id]
    
    df <- df[order(df$Gn, df$Tx, df$Pos, df$Type), ]

    rownames(df) <- NULL
    df[["order"]] <- 1:nrow(df)
    df
    
}

## creates a data frame containing all edges when a data frame of
## vertices is provided
edge <- function(v) {
    
    ## site type
    type <- v$Type
    vid <- v$Vid
    tx <-  v$Tx
    oo <- v$order
    
    ## gene name
    gn <- factor(v$Gn, levels=unique(v$Gn))
    
    ## root indices
    idx0 <- type=="0"
    
    ## leave indices
    idx1 <- type=="1"
    
    ## *[ -> 0->[
    ## start index
    idx <- type == "["
    
    ## create start point (vertex) for the edges of each gene
    ## the from vetor contains al start vertizes indices of all genes
    ## retrive root vertex idices
    from <- uunlist(split(vid[idx0], gn[idx0]))
    of <- uunlist(split(oo[idx0], gn[idx0]))
    
    
    ## start site of the transcripts
    to <- split(vid[idx], gn[idx])
    ot <- split(oo[idx], gn[idx])
    
    ## from root to the individual transcript start sites for each
    ## gene
    df0 <- data.frame(From=rep(from, sapply(to, length)),
                      To=uunlist(to), of=rep(of, sapply(to, length)),
                      ot=unlist(ot), oo=oo[idx])
    
    ## ]*, -> ]->1
    ## index of transcript ends
    idxr <- type == "]"
    from <- split(vid[idxr], gn[idxr])
    to <- uunlist(split(vid[idx1], gn[idx1]))
    
    of <- split(oo[idxr], gn[idxr])
    ot <- uunlist(split(oo[idx1], gn[idx1]))
    
    ## from the transcript ends to leaf
    df1 <- data.frame(From=uunlist(from),
                      To=rep(to, sapply(from, length)), of=unlist(of),
                      ot=rep(ot, sapply(from, length)), oo=oo[idxr])
    ## all others
    
    ## vertices between root and leaves
    idx <- !(idx0 | idxr | idx1)
    
    ## extract all vids except the last one 
    from <- head(vid, -1)[idx]
    of <- head(oo, -1)[idx]
    
    ## extract all vids except the first one
    to <- tail(vid, -1)[idx]
    ot <- tail(oo, -1)[idx]
    
    
    ## combine the edges from root to the first node
    ## from the last node to the leaves
    ## and all other not including root or leave
    df <- rbind(data.frame(From=from,To=to, oo=oo[idx], of=of, ot=ot),
                df0, df1)
    
    ## order edges
    df <- df[order(df$oo, df$From, df$To),]
    rownames(df) <- NULL
    
    ## create an edge id
    df$Eid <- seq_len(nrow(df))
    
    ## check for duplicates of edges
    d <- duplicated(df[,c("From", "To")])
    
    ## remove duplicates
    df[!d, c("Eid", "From", "To", "of", "ot")]
}









