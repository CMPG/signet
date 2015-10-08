#' Build adjacency matrix of a graph
#'
#' This function generates the adjacency matrix from a \verb{graphNEL} object.
#' @param pathway A \verb{graphNEL} object.
#' @param directed If \verb{TRUE}, edges direction will be considered.
#' @keywords adjacency matrix, graph
#' @return A n x n matrix, n is the number of nodes in the graph.
#' @examples
#' data(keggPathways)
#' adj<-getAdjacencyMatrix(keggPathways[[1]])

getAdjacencyMatrix<-function(pathway,
                             directed = FALSE)
{
  require(graph)
  x<-edges(pathway)
  GList <- names(x)
  adjMatrix <- matrix(0,length(GList),length(GList),dimnames=list(GList,GList))

  if(!directed)
  {
    for(i in GList)
    {
      adjMatrix[i,x[[i]]]<-1
      adjMatrix[x[[i]],i]<-1
    }
  }
  else
  {
    for(i in GList)
    {
      adjMatrix[i,x[[i]]]<-1
    }
  }

  rownames(adjMatrix)<-colnames(adjMatrix)<-GList

  return(adjMatrix)
}
