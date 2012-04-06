\name{spliceGraph}

\alias{spliceGraph}
\alias{spliceGraph,TranscriptDb-method}


\title{
  Creates a splicing graph from a \link{TranscriptDb} object
}

\description{
  Extracts the edges and the corresponding types of splicing events from a splicing graph 
  structure created based on the provided \link{TranscriptDb} object.
}

\usage{
spliceGraph(txdb, ...)
\S4method{spliceGraph}{TranscriptDb}(txdb, getEvents=TRUE, genes=getGH(txdb), ...)
}

\arguments{  
  \item{txdb}{
     A \link{TranscriptDb} object.
  }

  \item{getEvents}{
     Logical, if TRUE splice event information is provided. If the users interest 
     lies only on the splice graph edges without considering the types of splicing 
     events this argument should be set to FALSE to shorten the computation time.
  }
  
  \item{genes}{
    A character vector representing the ENTREZ gene ids used for creating the splicing 
    graph. If no genes are specified all genes contained in the \link{TranscriptDb} 
    object are used.
  }
  \item{...}{
    Arguments to be passed to or from methods.
  }
}

\details{
  This is the main function for creating and processing splicing graph structures created 
  from the provided \link{TranscriptDb} object. 
  Splicing graphs in general are efficient representations 
  of annotated transcript variants of genes and can provide also information about the 
  underlying splicing events generating those variants.
  The main elements of splicing graphs are called vertices and edges. 
  Vertices have fixed fixed genomic coordinates and represent basically potential splice sites.
  Edges mimic the corresponding intronic and exonic structure of the transcripts. 
  Since not all potential splice sites are 
  alternative splice sites splice graphs can be simplified by collapsing edges and vertices 
  not associated with alternative splicing events.
  The collapsed representation of the splice graph may contain bubbles between two
  remaining vertices. Such bubbles only occur if there exist more than one structurally different 
  transcript parts flanked by the two vertices, resulting from different choices 
  of the splicing machinery.
  The number of structurally different parts, and their properties represented by the
  edges and vertices in the non-collapsed graph allow conclusions about the 
  underlying splicing mechanism leading to the individual observed transcript variants. 
  Finally the function provides a \link{GRangesList} object containing the collapsed 
  edges of the splice graph and their associated exons. 
  The \link{metadata} slot provides also information about the type of splice 
  events occurring in the individual genes and their associated bubble IDs.
  It is important to note that the exons per edges in the 
  \link{GRangesList} are not the original exons defined within the \link{TranscriptDb} 
  object. The original exons, the genes and the new exons can be retrieved from the
  resulting \link{GRangesList} object. See the example code below.
}

\value{ 
  Returns a \link{GRangesList} object containing the collapsed edges of the splicing
  graph. This object contains the new disjoint exons, the gene ids, the original
  exon ids and the observed splicing events if the getEvents argument is set TRUE. 
}

\author{
  M. Carlson, M. Morgan, D. Bindreither
}

\references{
  Heber, S., Alekseyev, M., Sze, S., Tang, H., and Pevzner, P. A. 
  \emph{Splicing graphs and EST assembly problem}
  Bioinformatics
  Date: Jul 2002
  Vol: 18
  Pages: S181-S188

  Sammeth, M. (2009) 
  \emph{Complete alternative splicing events are bubbles in splicing graphs}
  J. Comput. Biol.
  Date: Aug 2009
  Vol: 16
  Pages: 1117-1140
}

\seealso{
  \link{TranscriptDb}, \link{GRangesList}
}

\examples{
  if(interactive()) {
    library(TxDb.Hsapiens.UCSC.hg18.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene

    ## creates the splice graph and retrieves the exons
    ## associated with the individual edges
    gn <- c("100", "100033417", "3726")

    sgModel <- spliceGraph(txdb, genes=gn)
    sgModel
    
    ## get the new exon ids
    values(unlist(sgModel))[["disJ_exon_id"]]
    
    ## get the gene ids
    values(sgModel)[["gene_id"]]
    
    ## get the original exon ids
    values(unlist(sgModel))[["exon_ids"]]

    ## access the splice events via the metadata slot
    spliceEvents <- metadata(sgModel)$spliceEvents

    head(spliceEvents)
    
    ## quantify the splice events
    sp <- split(spliceEvents$bubbleCode, 
                factor(spliceEvents$eventType))

    eventSum <- data.frame(event=names(sp), Nr=elementLengths(sp),
                           row.names=NULL)

    head(eventSum[order(eventSum$Nr, decreasing=TRUE), ])

    ## assignment of edges to bubbles
    values(sgModel)

  }
}
