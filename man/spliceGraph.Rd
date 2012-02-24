\name{spliceGraph}

\alias{spliceGraph}
\alias{spliceGraph,TranscriptDb-method}


\title{
  Create a splicing graph from an \link{TranscriptDb} object
}

\description{
  Extracts the edges of a splicing graph structure which is generated based
  on a given \link{TranscriptDb} object.
}

\usage{
spliceGraph(txdb, ...)
\S4method{spliceGraph}{TranscriptDb}(txdb, genes=getGH(txdb), ...)
}

\arguments{  
  \item{txdb}{
     A \link{TranscriptDb} object.
  }
  \item{genes}{
    A character vector representing the ENTREZ gene ids used for creating the splicing graph.
    If no genes are specified all genes in the \link{TranscriptDb} object are used.
  }
  \item{...}{
    Arguments to be passed to or from methods.
  }
}

\details{
  This is the main function for creating a splicing graph structure based on a provided 
  \link{TranscriptDb} object. Splicing graphs in general are efficient representations of 
  gene structures and associated alternative splicing information. The main elements of a 
  splicing graph are called vertices and edges. Vertices (or nodes) basically represent 
  potential splicing events and the edges mimic the corresponding 
  introns and exons of all known transcipt variants. Since not all splice sites are 
  alternative splice sites the graph can be simplified by collapsing edges and the vertices 
  where no evidence for alternative splicing is present. Finally the function provides a 
  \link{GRangesList} object containing the collapsed edges and their associated exons.
  Each edge contains at least one exon. Keep in mind that the exons per edges in the 
  \link{GRangesList} are not the original exons defined within the \link{TranscriptDb} 
  object. The original exons, the genes and the new exons can be rerieved from the
  resulting \link{GRangesList} object. See the example code below.
}

\value{ 
  Returns a \link{GRangesList} object containing the collapsed edges of the splicing
  graph. This object contains the new disjoint exons, the gene ids as well as the original
  exon ids. 
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

    ## creates the splice graph structure and retrieves the exons
    ## associated with the individual edges
    exonsByEdges <- spliceGraph(txdb, genes=c("1", "10"))
    exonsByEdges
    
    ## get the new exon ids
    values(unlist(exonsByEdges))[["disJ_exon_id"]]
    
    ## get the gene ids
    values(exonsByEdges)[["gene_id"]]
    
    ## get the original exon ids
    values(unlist(exonsByEdges))[["exon_ids"]]
  }
}
