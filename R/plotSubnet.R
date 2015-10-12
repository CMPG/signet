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
#' plotSubnet()

plotSubnet<-function(path,liste)
{
  subs<-as.character(liste)

  col<-c(rep("red",length(subs)))
  nAttrs<-list()
  nAttrs$fillcolor <- col
  nAttrs$height <- nAttrs$width <- rep("0.4", length(graph::nodes(path)))
  names(nAttrs$width)<-names(nAttrs$height)<-graph::nodes(path)
  names(nAttrs$fillcolor)<-c(subs)

  return(plot(path, y="neato", nodeAttrs = nAttrs))
}
