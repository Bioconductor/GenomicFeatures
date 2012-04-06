uunlist <- GenomicFeatures:::uunlist

test_fastDisjoinOnGRangesList <- function(){
    G1 <- GRanges(seqnames="chrX",
                  ranges=IRanges(start=c(1), end=c(10)),
                  strand="+")
    G2 <- GRanges(seqnames="chrX",
                  ranges=IRanges(start=c(6), end=c(12)),
                  strand="+")
    GR.test <- GRangesList(G1=G1, G2=G2)
    GR.new <- GenomicFeatures:::fastDisjoinOnGRangesList(GR.test)
    checkIdentical(elementLengths(GR.test), elementLengths(GR.new))
}

test_vertex <- function(){
    ## Test case 2: Disjoining often leads to exons which consist only of one NT.
    ## Because this could cause problems in generating the splicing graph since acceptor,
    ## donor start and end type can fall into on single vertices.
    ## T1 ||||||......
    ## T2 ......|.....
    ## T3 .......|||||

    T1 <- GRanges("chrX", IRanges(1,6), strand="+")
    values(T1)[["exon_id"]] <- 1111L
    values(T1)[["exon_rank"]] <- 2L
    T2 <- GRanges("chrX", IRanges(7,7), strand="+")
    values(T2)[["exon_id"]] <- 2222L
    values(T2)[["exon_rank"]] <- 1L
    T3 <- GRanges("chrX", IRanges(8,12), strand="+")
    values(T3)[["exon_id"]] <- 3333L
    values(T3)[["exon_rank"]] <- 1L
    
    test <- GRangesList( T1, T2, T3)
    names(test) <- c(1,2,3)
    
    test1 <- uunlist(GRangesList(GRanges("chrX", IRanges(1,6), strand="+"),
                                 GRanges("chrX", IRanges(7,7), strand="+"),
                                 GRanges("chrX", IRanges(8,12), strand="+")))
    
    values(test1)[["tx_id"]] <- c(1,2,3)
    seqlengths(test1) <- 12

    exsByTxs.clean <- test
    seqlengths(exsByTxs.clean) <- 12
    txsByGene.clean <- GRangesList(test1)
    names(txsByGene.clean) <- 1
    seqlengths(txsByGene.clean) <- 12


    txIdsInGeneOrder <- values(uunlist(txsByGene.clean))[["tx_id"]]

    ## Retrieve the gene ids from the transcripts by genes object
    gnIdsTxsByGns <-
      rep(names(txsByGene.clean), elementLengths(txsByGene.clean))

    ## Create a map of the gene ids (elements) to transcript ids
    ## (element names)
    geneTxMap <- setNames(gnIdsTxsByGns,
                          values(uunlist(txsByGene.clean))[["tx_id"]])

    exsByTxs.clean <- exsByTxs.clean[as.character(txIdsInGeneOrder)]
    
    dJExsByGns <-
      GenomicFeatures:::disjoinExonsByGene(exsByTxs.clean,
                                           txsByGene.clean)

    dJExsByTxs <- GenomicFeatures:::disjoinExtx(exsByTxs.clean,
                                                dJExsByGns, geneTxMap)

    siteMap <- setNames(as.raw(0:6),
                        c(" ", "0", "[", "^", "-", "]", "1"))
    bm <- any(duplicated(GenomicFeatures:::vertex(dJExsByTxs,
                                                  dJExsByGns,
                                                  siteMap)$Vid))
    
    checkIdentical(FALSE, bm)
       
}


test_bubble_processing <- function() {

    ## Test case 3: Create an example bubble

    ##     o - - - - - - o
    ##    -    -       -    
    ##  -        -   -
    ## o  - - - - - o 

  
    ## create an artifical example of 3 intersecting bubbles
    Ex1 <- GRanges("chrX", IRanges(1,6), strand="+")
    Ex2 <- GRanges("chrX", IRanges(8,12), strand="+")
    Ex3 <- GRanges("chrX", IRanges(16,20), strand="+")
    Ex4 <- GRanges("chrX", IRanges(24,26), strand="+")
    Ex5 <- GRanges("chrX", IRanges(30,34), strand="+")

    Tx1 <- c(Ex1, Ex4)
    values(Tx1)[["exon_id"]] <- c(1,4)
    values(Tx1)[["exon_rank"]] <- c(1,4)
    Tx2 <- c(Ex1, Ex3, Ex5)
    values(Tx2)[["exon_id"]] <- c(1,3,5)
    values(Tx2)[["exon_rank"]] <- c(1,3,5)
    Tx3 <- c(Ex2, Ex5)
    values(Tx3)[["exon_id"]] <-  c(2,5)
    values(Tx3)[["exon_rank"]] <- c(2,5)
    
    exsByTxs.clean <- GRangesList(Tx1, Tx2, Tx3)
    seqlengths(exsByTxs.clean) <- 34
    names(exsByTxs.clean) <- 1:3


    Tx1 <- GRanges("chrX", IRanges(min(start(Tx1)),
                                   max(end(Tx1))), strand="+")
    Tx2 <- GRanges("chrX", IRanges(min(start(Tx2)),
                                   max(end(Tx2))), strand="+")
    Tx3 <- GRanges("chrX", IRanges(min(start(Tx3)),
                                   max(end(Tx3))), strand="+")

    Tx.all <- c(Tx1, Tx2, Tx3)
    values(Tx.all)[["tx_id"]] <- 1:3
    txsByGene.clean <- GRangesList(Tx.all)
    seqlengths(txsByGene.clean) <- 34
    names(txsByGene.clean) <- "Gn1"


    ## construct the splice graph
    dJExsByGns <-
      GenomicFeatures:::disjoinExonsByGene(exsByTxs.clean,
                                           txsByGene.clean)
    
    geneTxMap <- setNames(rep("Gn1", 3), 1:3)    
    dJExsByTxs <-
      GenomicFeatures:::disjoinExtx(exsByTxs.clean, dJExsByGns,
                                    geneTxMap)
   
    siteMap <- setNames(as.raw(0:6), c(" ", "0", "[", "^",
                                       "-", "]", "1"))
   
    vertices <-
      GenomicFeatures:::vertex(dJExsByTxs, txsByGene.clean, siteMap)

    vertices.orig <- GenomicFeatures:::vertex(exsByTxs.clean,
                            txsByGene.clean, siteMap, exID="exon_id")
    code.orig <- apply(vertices.orig[, c("Pos", "Gn", "Tx")], 1,
                       paste, collapse=":")

    code.dj <- apply(vertices[, c("Pos", "Gn", "Tx")], 1,
                     paste, collapse=":")

    vertices$origType <- ifelse(code.dj %in% code.orig,
                                as.character(vertices$Type),
                                paste(vertices$Type, ":", sep=""))

    edgesPrelim <- GenomicFeatures:::edge(vertices)
    
    edges <- GenomicFeatures:::rmDupEdges(edgesPrelim)
    
    etscript <- GenomicFeatures:::getEtscript(edgesPrelim)

    verticesIod <- GenomicFeatures:::iodeg(vertices, edges[, seq(1, 5)])

    cEdges <- GenomicFeatures:::collapseEdge(verticesIod, edges)
    cVertices <- GenomicFeatures:::collapseVertex(verticesIod)
    vtscript <- GenomicFeatures:::vertexTranscripts(vertices)
    cVtscript <- vtscript[as.character(cVertices$Vid)]
    
    cEtscript <-
      GenomicFeatures:::edgeTranscripts(vtscript, edges, cEdges)

    ## test bubble extraction
    allBub <- GenomicFeatures:::bubbles(cEtscript, cVertices, cEdges,
                                        cVtscript, vertices,
                                        origVids=FALSE )

    ## test bubble to edge mapping
    l <- GenomicFeatures:::collapseEdgeX(verticesIod, edges)
    temp <- GenomicFeatures:::edgeToBubble(allBub, l)
    bub <- any(!names(temp) %in% cEdges$Eid)
    temp <- GenomicFeatures:::edgeToBubblePart(allBub, l)
    bubp <- any(!names(temp) %in% cEdges$Eid)

    ## check if bubbles could be extracted and if the collapsed edges
    ## to vertices map is correct
    identical(FALSE, any(bubp, bub, is.null(allBub)))
}
