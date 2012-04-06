setMethod("spliceGraph", signature(txdb="TranscriptDb"),
    function(txdb, getEvents=TRUE, genes=getGH(txdb), ...)
    {
        ## Get transcripts per gene
        txsByGene <- transcriptsBy(txdb, "gene")
        
        ## Subset to a small number of genes
        if(class(genes) == "character" & length(genes) > 0) {
          if(any(! as.character(genes) %in% names(txsByGene)) &
             any(! getGH(txdb) %in% genes))
            stop("genes have to match the gene ids in txdb!")
          inclEntrezGeneIds <- as.character(genes)
          txsByGene <- txsByGene[names(txsByGene) %in% inclEntrezGeneIds]
        }
        
        ## Clean transcripts by gene from strange transcripts and from
        ## transcripts arising from strange chromosomes
        txsByGene <- txclean(txsByGene)
        
        
        exsByTxs <- exonsBy(txdb, "tx")
        rm(txdb)
        txIds <- values(unlist(txsByGene,
                               use.names=FALSE))[["tx_id"]]
        
        exsByTxs <- exsByTxs[as.character(txIds)]
        
        ## Clean up the exons  by transcripts
        exsByTxs <- txclean(exsByTxs)
        exIdsTxOrder <- values(uunlist(exsByTxs))[["exon_id"]]
        txIdsExsByTxs <- rep(names(exsByTxs),
                             elementLengths(exsByTxs))
        
        
        ## Combine the exon ids and transcript ids to a generate a 
        ## look up data frame which expanded by adding the associated
        ## gene ids later
        exToTx <- data.frame(Tx=txIdsExsByTxs, Ex=exIdsTxOrder)
        
        ## Generating a gene to transcript mapping
        ## Get all transcript ids of the transcripts by gene object
        txIdsInGeneOrder <-
          values(uunlist(txsByGene))[["tx_id"]]
        
        ## Create a map of the gene ids (elements) to transcript
        ## ids (element names)
        gnIdsTxsByGns <- rep(names(txsByGene),
                             elementLengths(txsByGene))
        
        geneTxMap <- setNames(gnIdsTxsByGns, txIdsInGeneOrder)
        
        ## Include gene ids into the look up data frame exons to
        ## transcripts
        
        exToTxToGn <- exToTx
        exToTxToGn$Gn <- geneTxMap[exToTx$Tx]
        rm(exToTx)
        
        ## ensure that the ordering of the exons per transcript is
        ## congruent with the order in the transcripts by genes object
        exsByTxs <-
          exsByTxs[as.character(txIdsInGeneOrder)]
        
        ## disjoin exons per gene to remove congruent overlaps
        ## keep in mind that disjoining discards exon ids coming from 
        ## the TxDb object and creates new exon ids
        dJExsByGns <-
          disjoinExonsByGene(exsByTxs, txsByGene)
        
        
        ## disjoin exons per transcript, here was the first bug caused
        ## by not taking the gene ID into acount during disjoining
        ## which caused exon mix-ups between genes who shared common 
        ## exons identical regarding there genomic coordinates
        
        dJExsByTxs <-
          disjoinExtx(exsByTxs, dJExsByGns, geneTxMap)
        
        ### Create the splicing graph
        
        ## Mapping of the splice site types
        siteMap <- setNames(as.raw(0:6), c(" ", "0", "[", "^",
                                           "-", "]", "1"))
        
        txId <- values(uunlist(txsByGene))[["tx_id"]]
        
        ## Create the vertices of the graph (data.frame)
        vertices <- vertex(dJExsByTxs[as.character(txId)],
                           txsByGene, siteMap)


        ## create original vertices before dissjoining to
        ## carry forward the original gene model information
        vertices.orig <- vertex(exsByTxs, txsByGene,
                                siteMap, exID="exon_id")


        ## get original vertices to identify alternative
        ## acceptor and donor sites this is éssential to differentiate
        ## between exon skips and alternative donor sites later on
        code.orig <- paste(vertices.orig$Pos, vertices.orig$Gn,
                           vertices.orig$Tx, sep=":")

        code.dj <- paste(vertices$Pos, vertices$Gn, vertices$Tx,
                         sep=":")

        vertices$origType <- ifelse(code.dj %in% code.orig,
                                    as.character(vertices$Type),
                                    paste(vertices$Type, ":", sep=""))
        
        ## create the edges data frame of the splicing graph
        edgesPrelim <- edge(vertices)

        ## remove duplicated edges
        etscript <- getEtscript(edgesPrelim)
        edges <- rmDupEdges(edgesPrelim)
       
        rownames(vertices) <- vertices$order

        
        ## compute in and out degree of the individual vertices
        verticesIod <-
          iodeg(vertices, edges[, c("Eid", "From", "To", "of", "ot")])
        
        ## collapse the edges, to only keep those where at least one
        ## associated vertex has outdegree or indegree not equal to 1
        cEdges <- collapseEdge(verticesIod, edges)

        ## get transcripts associated with each vertex
        vtscript <- vertexTranscripts(vertices)

        ## transcripts associated with the collapsed edges
        cEtscript <- edgeTranscripts(vtscript, edges, cEdges)
        
        ## collapse the vertices, we keep only those which have
        ## outdegree or indegree not equal to 1
        cVertices <- collapseVertex(verticesIod)

        ## transcripts associated with the collapses vertices 
        cVtscript <- vtscript[as.character(cVertices$Vid)]
        
        ## retrieve all exons per edges
        eexons <- collapseEdgeExons(verticesIod, edges)
        
        ## retrieve all exons per edges of all edges
        exsByEdges <- exonsByEdge(dJExsByGns, eexons)
        
        ## create a edge to gene mapping for comparison purposes
        vid <- cEdges[match(names(exsByEdges), cEdges$Eid), "From"]
        gn <- cVertices[match(vid, cVertices$Vid),"Gn"]
        
        edgeToGene <- data.frame(Edge=names(exsByEdges), Gene=gn)
        rownames(edgeToGene) <- edgeToGene$Edge
        
        ## assign gene names to the individual edges
        values(exsByEdges)[["gene_id"]] <-
          as.character(edgeToGene[names(exsByEdges), ]$Gene)

        if(getEvents) {
            allBubbles <- bubbles(cEtscript, cVertices, cEdges,
                                  cVtscript, vertices,
                                  origVids=FALSE )
            
            sc <- spliceCodeCharacterList(allBubbles)
            
            strRes <- spliceCodeToString(sc)
            
            ## get the event description
            fn <- system.file("extdata", "events.Rda",
                                  package="GenomicFeatures")

            load(fn)

            l <- collapseEdgeX(verticesIod, edges)
            
            values(exsByEdges)[["bubble_ids"]] <-
              CharacterList(edgeToBubble(allBubbles, l)[names(exsByEdges)])

            values(exsByEdges)[["bubblePart_ids"]] <-
              CharacterList(edgeToBubblePart(allBubbles, l)[names(exsByEdges)])

            df <- DataFrame(bubbleCode=names(sc),
                            code=sc,
                            eventType=unname(events[strRes]),
                            gn=elementMetadata(sc)$Gn,
                            bubbleNr=as.integer(elementMetadata(sc)$bubble)
                            )
       
            metadata(exsByEdges) <- list(spliceEvents=df)
        }
        
        exsByEdges    
  }
)
