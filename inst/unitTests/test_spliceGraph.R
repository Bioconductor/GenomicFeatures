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
    gnIdsTxsByGns <- rep(names(txsByGene.clean), elementLengths(txsByGene.clean))

    ## Create a map of the gene ids (elements) to transcript ids (element names)
    geneTxMap <- setNames(gnIdsTxsByGns, values(uunlist(txsByGene.clean))[["tx_id"]])

    exsByTxs.clean <- exsByTxs.clean[as.character(txIdsInGeneOrder)]
    
    dJExsByGns <- GenomicFeatures:::disjoinExonsByGene(exsByTxs.clean, txsByGene.clean)

    dJExsByTxs <- GenomicFeatures:::disjoinExtx(exsByTxs.clean, dJExsByGns, geneTxMap)

    siteMap <- setNames(as.raw(0:6), c(" ", "0", "[", "^", "-", "]", "1"))
    
    checkIdentical(FALSE, any(duplicated(GenomicFeatures:::vertex(dJExsByTxs, dJExsByGns, siteMap)$Vid)))
       
}

