#' plotSubnet
#'
#' This function allows you to get the adjacency matrix from a graph object.
#' @param pathway A graph object.
#' @param liste Set to TRUE if you want the graph to be directed.
#' @keywords adjacency matrix
#' @export
#' @examples
#' #get data
#' require(signet)
#' #plotSubnet()

plotSubnet<-function(pathway,liste)
{
  subs<-as.character(liste)

  col<-c(rep("red",length(subs)))
  nAttrs<-list()
  nAttrs$fillcolor <- col
  nAttrs$height <- nAttrs$width <- rep("0.4", length(graph::nodes(pathway)))
  names(nAttrs$width)<-names(nAttrs$height)<-graph::nodes(pathway)
  names(nAttrs$fillcolor)<-c(subs)

  return(plot(pathway, y="neato", nodeAttrs = nAttrs))
}
