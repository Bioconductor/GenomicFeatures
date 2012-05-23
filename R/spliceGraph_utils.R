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
                                                                      
txclean <- function(tx, seqlevels=paste("chr", seq(1,21), sep=""), ...){

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
    ## Add in-degree and out-degree
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
    tscript <- lapply(split(v$Tx[idx], f), unique)
    
    ## Give the list of transcripts per genes as name the root and leaf
    ## vertices idex. Assign to the individual vertices the
    ## corresponding transcript ids
    test <- c(setNames(tscript, v$Vid[v$Type=="0"]),
              setNames(tscript, v$Vid[v$Type=="1"]),
              split(v$Tx[idx], v$Vid[idx]))
    
    test
}

edgeTranscripts <- function(vtscript, e, e2) {

    ue <- e[e$From %in% e2$From,]
  
    ## extract for each remaining edge those transcripts which
    ## tx id overlap of "from" node and "to" node of the edges
    etscript <- Map(intersect, vtscript[as.character(ue$From)],
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

    ## remove duplications
    nn <- rep(names(orig), elementLengths(orig))
    orig.flat <- uunlist(orig)
    idx <- !duplicated(paste(nn, orig.flat, sep=":"))
    nn <- nn[idx]
    orig.flat <- orig.flat[idx]
    orig <- split(orig.flat, nn)
    
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
    ot <- e$ot
    rownames(v) <- NULL

    end <- unique(v$Vid[v$Indeg!=1 | v$Outdeg != 1])
    
    start <- unique(v$Vid[v$Outdeg==1])

    ## same as collpase edges
    x <- start[!start %in% end]
    
    while (length(x)) {
        i <- match(x, to)
        j <- match(x, e$From)
        to[i] <- x <- to[j]
        edge[j] <- edge[i]
        x <- x[!x %in% end]
    }

    ## to conatins the new ends of the exons
    idx <- match(e$From, v$Vid)

    ## ok until here
    v$Type[idx]
    
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
vertex <- function(extx, txgn, map, exID="disJ_exon_id") {
    
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
    eex <- as.character(values(ex)[[exID]])

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

    ex.names <- values(ex)[[exID]]
    uni.idx <- ! duplicated(ex.names)    
    ex.uni <- ex[uni.idx]
    
    ex.widths <- width(ranges(ex.uni))
    names(ex.widths) <- values(ex.uni)[[exID]]
    
    ## create a exon data frame
    exdf <- data.frame(Pos=epos, Gn=tgn[etx], Tx=etx,
                       Ex=eex, Type=etype, stringsAsFactors=FALSE)

    ## check for one NT long exons, if one vertices is assiciated with
    ## such an exon keep the vertices anyway, very crucial part
    leng <- ex.widths[exdf$Ex]
    z <- paste(leng, exdf$Type, exdf$Pos, sep=":")
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
    return(df)
    
}

## creates a data frame containing all edges when a data frame of
## vertices is provided
edge <- function(v) {

    ## avoids duplication in the resulting outcome
    ## very crucial step in the code
    s <- ifelse(v$Type %in% c("[", "-"), "q", "r")
    s <- ifelse(v$Type %in% c("]", ""), "e", s)

    v <- v[v$keep == "K" | ! duplicated(paste(v$Gn, v$keep=="K",
             v$Pos, s, sep=":")), ]
    
    ## site type
    type <- v$Type
    vid <- v$Vid
    tx <-  v$Tx
    oo <- v$order
    oGn <- v$Gn
    
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
    genes <- uunlist(split(as.character(gn)[idx], as.character(gn)[idx]))

    
    ## from root to the individual transcript start sites for each
    ## gene
    df0 <- data.frame(From=rep(from, sapply(to, length)),
                      To=uunlist(to), of=rep(of, sapply(to, length)),
                      ot=unlist(ot), oo=oo[idx], tx="0", gn = genes)
    
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
                      ot=rep(ot, sapply(from, length)), oo=oo[idxr],
                      tx="N", gn=oGn[idxr])
    ## all others
    
    ## vertices between root and leaves
    idx <- !(idx0 | idxr | idx1)
    
    ## extract all vids except the last one 
    from <- head(vid, -1)[idx]
    of <- head(oo, -1)[idx]
    
    ## extract all vids except the first one
    to <- tail(vid, -1)[idx]
    ot <- tail(oo, -1)[idx]
    genes <- tail(gn, -1)[idx]

    txf <- tail(tx, -1)[idx]
    
    
    ## combine the edges from root to the first node
    ## from the last node to the leaves
    ## and all other not including root or leave
    df <- rbind(data.frame(From=from,To=to, oo=oo[idx], of=of, ot=ot,
                           tx=txf, gn=genes), df0, df1)
    
    ## order edges
    df <- df[order(df$oo, df$From, df$To),]
    rownames(df) <- NULL
    
    ## create an edge id
    df$Eid <- seq_len(nrow(df))

    df
}


rmDupEdges <- function(df) {
    
    ## check for duplicates of edges
    d <- duplicated(df[,c("From", "To", "gn")])
    
    df[!d, c("Eid", "From", "To", "of", "ot", "gn"), ]
  }

getEtscript <- function(df) {

    crit <- paste(df$From, df$To, df$gn, sep=":")

    edgeT <- split(as.character(df$tx), crit)

    d <- duplicated(df[,c("From", "To", "gn")])

    df <- df[!d, ]
    rownames(df) <- crit[!d]
    
    ## remove duplicates
    names(edgeT) <- df[names(edgeT),]$Eid

    ## check for root an leaf duplicates
    edgeT
}


## assigns the edge ids to the bubbles
edgeToBubble <- function(b, vidToEdge) {

    ## create a map of original v to edges
    v <- split(rep(names(vidToEdge), elementLengths(vidToEdge)),
          as.character(uunlist(vidToEdge)))

    ## look up the edges associated with each v
    v <- v[names(v) %in% unique(as.character(b$origVid))]        
    v <- v[match(as.character(b$origVid), names(v))]

    ## create the bubble code
    e <- rep(paste(b$Gn, b$bubble, sep=":"), elementLengths(v))

     ## unify the bubble to edge mapping
    v <- as.character(uunlist(v))
    uni <- ! duplicated(paste(e, v, sep=":"))
    split(e[uni], v[uni])
}

collapseBubbleParts <- function(eToBlp) {
    ## collapse the bubble parts for the edges if they
    ## are redundant
    bub <- uunlist(eToBlp)
    en <- rep(names(eToBlp), elementLengths(eToBlp))
    
    
    temp <-
      t(simplify2array(strsplit(bub, ":")))[, seq(1,2)]
    temp <- paste(temp[,1], temp[,2], sep=":")
    
    
    d <- data.frame(cbind(en, bub, temp,
                          be=paste(temp, en, sep="-")))
    d <- d[order(d$temp, d$en, d$bub), ]
    
    idx <- match(unique(d$be), d$be)
    
    bpmap <- setNames(as.character(d$bub[idx]),
                      as.character(d$be[idx]))

    
   

    map <- data.frame(nn = as.character(bpmap[as.character(d$be)]),
                      oo = as.character(d$bub))

    map <- map[as.character(map$nn) != as.character(map$oo),
               , drop=FALSE]

    map <- setNames(as.character(map$nn), as.character(map$oo))
    
    

    repl <- as.character(map[as.character(d$bub)])
    
    d$collapse <- ifelse(!is.na(repl), repl, as.character(d$bub) )
    
    d <- d[! duplicated(paste(d$en, d$collapse, sep=".")), ]
    
    split(d$collapse, d$en)
}


## assigns the edge ids to the bubble parts
edgeToBubblePart <- function(b, vidToEdge) {

    ## create a map of original v to edges
    v <- split(rep(names(vidToEdge), elementLengths(vidToEdge)),
          as.character(uunlist(vidToEdge)))

    ## look up the edges associated with each v
    v <- v[names(v) %in% unique(as.character(b$origVid))]
    v <- v[match(as.character(b$origVid), names(v))]

    ## create the bubble part code
    e <- rep(paste(b$Gn, b$bubble, b$TxN, sep=":"),
             elementLengths(v))

    ## unify the bubble part to edge mapping
    v <- as.character(uunlist(v))
    uni <- ! duplicated(paste(e, v, sep=":"))
    collapseBubbleParts(split(e[uni], v[uni]))
}


## creates a compressed CharacterList of bubble splice codes
createSpliceCode <- function(codeV) {
    
    ## create ids for odering purposs
    codeV$bubbleTxID <- paste(codeV$Gn, codeV$bubble, codeV$Tx,
                              sep=":")
    codeV$bubbleID <- paste(codeV$Gn, codeV$bubble, sep=":")
    
    txGrp <- match(codeV$bubbleTxID, codeV$bubbleTxID)
    codeV$TxN <- codeV$Vid[txGrp]
    times <- elementLengths(split(txGrp, txGrp))
    
    ## ensure the right ordering of the vids
    ens <- c(diff(txGrp), 1)
    idx.last <- which(ens != 0)
    lastTxVid <- codeV$Vid[rep(idx.last, times)]
    
    ## order according genes, bubbles, tx start and tx end
    codeV <- codeV[order(paste(codeV$Gn, codeV$bubble, sep=":"),
                         codeV$TxN, lastTxVid), ]
    
    ## numerate artificial Tx accordin to the ordering employed above
    d <- c(0,diff(as.integer(codeV$Tx)))
    codeV$TxN <- cumsum(ifelse(d==0, 0, 1)) + 1
    
    bubGrp <- match(codeV$bubbleID, codeV$bubbleID)
    codeV$TxN <- codeV$TxN - codeV$TxN[bubGrp] + 1
    codeV
}

## maps the original vertices ids to the collapsed edges
collapseEdgeX <- function(v, e){
    ## collapse edge with Indeg == Outdeg == 1 
    to <- e$To

    edge <- e$Eid
    
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
        edge[j] <- edge[i]
        x <- x[!x %in% end]
    }
    
    i <- e$From %in% c(v$Vid[v$Outdeg!=1], v$Vid[v$Indeg!=1])

    split(e$From[!i], edge[!i])

}

## creates a compressed CharacterList of bubble splice codes
createSpliceCode <- function(codeV) {
    
    ## create ids for odering purposs
    codeV$bubbleTxID <- paste(codeV$Gn, codeV$bubble, codeV$Tx,
                              sep=":")
    codeV$bubbleID <- paste(codeV$Gn, codeV$bubble, sep=":")
    
    txGrp <- match(codeV$bubbleTxID, codeV$bubbleTxID)
    codeV$TxN <- codeV$Vid[txGrp]
    times <- elementLengths(split(txGrp, txGrp))
    
    ## ensure the right ordering of the vids
    ens <- c(diff(txGrp), 1)
    idx.last <- which(ens != 0)
    lastTxVid <- codeV$Vid[rep(idx.last, times)]
    
    ## order according genes, bubbles, tx start and tx end
    codeV <- codeV[order(paste(codeV$Gn, codeV$bubble, sep=":"),
                         codeV$TxN, lastTxVid), ]
    
    ## numerate artificial Tx accordin to the ordering employed above
    d <- c(0,diff(as.integer(codeV$Tx)))
    codeV$TxN <- cumsum(ifelse(d==0, 0, 1)) + 1
    
    bubGrp <- match(codeV$bubbleID, codeV$bubbleID)
    codeV$TxN <- codeV$TxN - codeV$TxN[bubGrp] + 1
    codeV
}

## removes duplicated edge information
rmNonInt <- function(codeV) {
    splitOption <- paste(codeV$Gn, codeV$bubble, codeV$Tx, sep=":")
    value <- paste(codeV$Gn, codeV$bubble, codeV$Vid, codeV$origType,
                   sep=":")
    
    e <-
      unlist(lapply(split(value, splitOption),
                    paste, collapse=" "))
    
    codeV[splitOption %in% names(e)[!duplicated(e)], ]  
}

getVertBetwST <- function(s, t, v, txsB, gn) {
    test <- data.frame(sources=match(as.integer(s[gn])+1, v$Vid),
                       sinks=match(t[gn], v$Vid))
    
    ## extract all v between source and sink without source and sink
    sel <-apply(test, 1, function(x) seq((x[1]), (x[2]-1)))
    v <- v[uunlist(as.vector(sel)), , drop=FALSE]
    
    ## remove vertices belonging to the source
    ## or sink 
    v <- v[!v$Vid %in% c(s,t), , drop=FALSE]
    
}

## adds dummy vids for txs with no present exon
addDummyVids <- function(v, txsB, txsBID) {
    noEx <- setNames(rep(names(txsB), elementLengths(txsB)), txsBID)
    bool <- which(! txsBID %in% v$Tx)
    
    if(length(bool) > 0) {
      noEx <- noEx[bool]
      dummy <- rep(NA, length(noEx))
      dummy0 <- rep(0, length(noEx))
      v <- rbind(v, data.frame(Pos=dummy, Gn=noEx,
                               Tx=names(noEx),
                               Type=dummy0,
                               Ex=dummy, keep=dummy,
                               Vid=dummy0, order=dummy,
                               origVid = dummy0,
                               origType=dummy0))
    }
    v
}

## group wise cumsum
fastListIdx <- function(gn, bool) {
    df <- data.frame(gn=as.character(gn), bool)
    cs <- cumsum(df$bool)
    aa <- match(df$gn, df$gn)
    cs2 <- c(0L, cs[-length(cs)])
    cs - cs2[aa]
}

## group wise sum
fastSum <- function(gn, bool) {
    cs <- fastListIdx(gn, bool)
    ro <- seq(length(cs), 1)
    rcs <- cs[ro]
    rgn <- gn[ro]
    unign <- unique(rgn)
    sort(unign[ rcs[match(unign, rgn)] > 0])
}

## extracts new bubbles
getBubbles <- function(txsB, v, t, s, gn, origV=FALSE, bc=1) {

    ## subset
    used <- intersect(names(s), names(t))
    t <- t[used]
    s <- s[used]
  
    ## get the transcript ids associated with the current bubbles
    txsBID <- uunlist(txsB)
   
    rownames(v) <- NULL
    
    v <- getVertBetwST(s, t, v, txsB, gn)

    # check for remaining transcripts
    txsB <- txsB[names(txsB) %in% unique(v$Gn)]
    txsBID <- uunlist(txsB)
    gn <- gn[gn %in% names(txsB)]
    v <- v[v$Gn %in% gn, ]

    ## if no vertices are remaining
    if(nrow(v) < 1) return(NULL)

    ## discard orig vertices ids
    if(! origV) {
      ## source vertices for each gene
      v <- v[order(v$Gn, v$Vid), ]
      idx <- as.integer(v[ match(gn, v$Gn), ]$Vid)
      names(idx) <- unique(v$Gn)
      v$origVid <- v$Vid
      v$Vid <- (as.integer(v$Vid) - idx[as.character(v$Gn)]) + 1
    }

    v <- v[order(v$Gn, v$Tx, v$Vid), ]
    
    ## if one of the txs has no exon at all
    ## introduce an artificial node per bubbble
    v <- addDummyVids(v, txsB, txsBID)

    ## set bubble number
    v$bubble <- bc

    ## remove duplicated information
    v <- rmNonInt(v[order(v$Gn, v$bubble, v$Tx, v$order), ])
    el <- elementLengths(lapply(split(v$Tx, v$Gn), unique))
    el <- el[el > 1]

    if(length(el) > 0) {

      v <- v[v$Gn %in% names(el), ]
      ## obtain the splice code
      ## Tx Vid VidType Tx Vid VidType .....
      q <- createSpliceCode(v[order(v$Gn, v$bubble, v$Tx, v$order), ])
      q[, c("Pos", "Gn", "Tx", "Vid", "order", "origType", "origVid",
            "bubble", "TxN")]

    } else return(NULL)
}

## represents the bubbles as character lists
spliceCodeCharacterList <- function(df.sc) {
    df.sc$bubbleID <- paste(df.sc$Gn, df.sc$bubble, sep=":")
  
    Tx <- split(df.sc$TxN, 1:nrow(df.sc))
    sign <- split(df.sc$origType, 1:nrow(df.sc))
    vids <- split(df.sc$Vid, 1:nrow(df.sc))
    
    perVidCode <- Map("c", Tx, vids, sign)
    names <- rep(df.sc$bubbleID, rep(3,nrow(df.sc)))
    codes <- split(uunlist(perVidCode), names)
    
    cl <- CharacterList(codes)
    m <- as(df.sc[match(unique(df.sc$bubbleID), df.sc$bubbleID),
               c("Tx", "Gn", "bubble")], "DataFrame")
    cl <- cl[paste(m$Gn, m$bubble, sep=":")]
    elementMetadata(cl) <- m
    cl
}

## creates tsring representation of the bubbles
spliceCodeToString <- function(cl){
    sapply(cl, paste, collapse=" ")
}


## retrieves the splicing code of all variants
## implements the algorithm proposed by sammeth 2009
bubbles <- function(cEtscript, cvt, cEdges, cVtscript, vertices,
                    d=2, origVids=FALSE) {

    ## reset rownames since the method is based on indexing
    rownames(cEdges) <- NULL
    rownames(cvt) <- NULL
    rownames(vertices) <- NULL
    
    vertices <- vertices[order(vertices$Gn, vertices$Vid), ]
    
    ## indentify sinks with indegree >= d
    t <- cvt$Indeg >= d
    t.idx <- fastListIdx(cvt$Gn, t)
    
    ## Assign the obtained sink indices to the vertices
    cvt <- data.frame(cvt, t.idx)
    
    ## Assign gene ID to edges data frame
    cEdges <- cEdges.orig <- data.frame(cEdges,
                         Gn = cvt$Gn[match(cEdges$To, cvt$Vid)])
    
    ## check for bubble less genes and remove them
    bubbleGns <- fastSum(cvt$Gn, t)

    ## initialize the data frame for storing the bubbles
    bubbles <- NULL
    
    # subset edges to only genes with bubbles
    cEdges <- cEdges[cEdges$Gn %in% bubbleGns, ]
    cvt <- cvt[cvt$Gn %in% bubbleGns, ]
    vertices <- vertices[vertices$Gn %in% bubbleGns, ]
    
    # initilize bubble counter
    bubble.cnt <- 1
    
    ## interate over sink vertices (t) in genomic order
    ## within each gene
    for(i in seq_len(max(t.idx))) {
      
      ## get first sinks for each gene
      tt <- match(unique(paste(i, cvt$Gn, sep=":")),
                  paste(cvt$t.idx, cvt$Gn, sep=":"))
    
      if(any(is.na(tt))) {
        tt <- tt[!is.na(tt)]
        cvt <- cvt[cvt$Gn %in% cvt$Gn[tt], ]
        tt <- match(unique(paste(i, cvt$Gn, sep=":")),
                    paste(cvt$t.idx, cvt$Gn, sep=":"))
    
        cEdges <- cEdges[cEdges$Gn %in% cvt$Gn, ]
      }
    
      ## Vids of the sinks and their associated genes
      tVids <- cvt[tt,]$Vid
      sinkGns <- cvt[tt,]$Gn
      names(tVids) <- sinkGns
      sinkGns <- rep(sinkGns,
                    elementLengths(cVtscript[as.character(tVids)]))
      
      ## get the edge ids of the incomming edges of the individual
      ## sinks
      inEdgeIDX <- cEdges$To %in% as.character(cvt[tt,]$Vid)
      sinkInEdges <- cEdges[inEdgeIDX, ]
    
      ## transcripts(t) per gene associated with the gene wise sinks
      temp <- cVtscript[as.character(cvt[tt,]$Vid)]
      chi <-
        split(uunlist(cVtscript[as.character(cvt[tt,]$Vid)]), sinkGns)
    
      ## initialize tuple collection
      C <- NULL
    
      ## keep the transcripts associated with the sinks
      sinkTxs <- chi

      ## order cEdges in genomic order regarding the from node
      cEdges <-
        cEdges[order(cEdges$From, cEdges$To, decreasing=FALSE), ]
    
      ## all vertices s < t with outdegreee >= 2
      times <- elementLengths(split(cvt$Gn, cvt$Gn))
      
      boolv2 <- cvt$Outdeg >= 2 & cvt$Vid < rep(cvt[tt,]$Vid, times)
    
      ## index the sources vertices within each gene
      s.idx <- fastListIdx(cvt$Gn, boolv2)
      cvt$s.idx <- s.idx
    
      ## initialize control over finished source genes
      cvt$S.fin <- FALSE
      
      # for all vertices s < t in reverse genomic order
      for(j in rev(seq_len(max(s.idx)))) {
        cvs <- cvt
        cEdgesS <- cEdges
        
        ## first source
        firstSs <-  match(unique(paste(j, cvs$Gn, sep=":")),
                          paste(cvs$s.idx, cvs$Gn, sep=":"))
    
        if(any(is.na(firstSs))) {
          firstSs <- firstSs[!is.na(firstSs)]
          cvs <- cvs[cvs$Gn %in% cvs$Gn[firstSs], ]
          firstSs <- match(unique(paste(j, cvs$Gn, sep=":")),
                           paste(cvs$s.idx, cvs$Gn, sep=":"))
          
          cEdgesS <- cEdgesS[cEdgesS$Gn %in% cvs$Gn, ]
        }
    
        
        # only active genes
        firstSs <- firstSs[!cvs$S.fin[firstSs]]
    
        ## Vids of the s and their associated genes
        sVids <- cvs[firstSs,]$Vid
        sourceGns <- cvs[firstSs,]$Gn
        names(sVids) <- sourceGns
        
        sourceGns <-
          rep(cvs[firstSs,]$Gn,
              elementLengths(cVtscript[as.character(sVids)]))
    
        sourceTxs <-
          split(unlist(cVtscript[as.character(cvs[firstSs,]$Vid)],
                       use.names=FALSE), sourceGns)
        
    
        # chiST = 0
        chiST <- as.list(rep(NA, length(names(chi))))
        names(chiST) <- names(chi)
    
        # outgoing edges
        edgidx <- cEdgesS$From %in% cvs[firstSs, ]$Vid
    
        # oder of the out going edges
        
        su.idx <- fastListIdx(cEdgesS$Gn, edgidx)
        cEdgesS$su.idx <- su.idx
    
        # control over source out edges
        cvs$su.fin <- FALSE
        
        # iterate over all s -> u edges
        for(k in seq_len(max(su.idx, na.rm=TRUE))) {
          firstSU <- match(unique(paste(k, cEdgesS$Gn, sep=":")),
                           paste(cEdgesS$su.idx, cEdgesS$Gn, sep=":"))
    
          firstSs <- firstSs[!cvs$su.fin[firstSs]]
    
          temp <- cEtscript[as.character(cEdgesS$Eid[firstSU])]
          gn <- as.character(rep(cEdgesS$Gn[firstSU],
                                 elementLengths(temp)))
          
          ## intersect Xst U X
          su <- split(uunlist(temp), gn)
          chiST.new <- Map(intersect, su , chi[names(su)])
    
          chiST[names(su)] <-
            Map(function(x,y) unique(c(x,y)[! is.na(c(x,y))]) ,
                chiST[names(su)], chiST.new[names(su)])
    
          ## remove chiST from chi per gene
          chi[names(su)] <- Map(function(x,y) y[!y %in% x],
                                chiST[names(su)], chi[names(su)])
        }
        
    
        ## Check for event extraction 
        extract <- elementLengths(chiST) >= d

        
        if(any(extract)) {
          outGn <- names(extract[extract])

          ## get bubble transcripts
          txsBubbles <- Map(intersect, sinkTxs[outGn],
                            sourceTxs[outGn])  

          ## check for new events
          txsBubbles <- txsBubbles[! txsBubbles %in% C]

          ## subset to genes where tx ar available
          outGn <- names(txsBubbles)
          
          if(length(txsBubbles) > 0) {
            newBubble <- getBubbles(txsBubbles, vertices, tVids,
                                    sVids, outGn, origVids, bubble.cnt)

            if(! is.null(newBubble)) {
              bubbles <- rbind(bubbles,newBubble)            
              bubble.cnt <- bubble.cnt + 1          
              C <- append(C, unname(txsBubbles))
            }
          }
        }
        
        ## X <- X U Xst
        chi <- Map("c", chi, chiST)
    
        ## check if transcripts(t) contained in transcripts(s) for
        ## all genes
        cond <- unlist(Map(function(x,y) all(x %in% y),
                           sinkTxs[names(sourceTxs)], sourceTxs ) )
        condPerV <- elementLengths(split(cvs$Gn,cvs$Gn))
          
        ## check for gene wise breaks
        cvs[cvs$Gn %in% names(cond), ]$S.fin <-
          rep(cond, condPerV[names(cond)])
        
        if(all(cond)) break
        
      }
    }
    bubbles
}
