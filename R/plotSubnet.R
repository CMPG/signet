#' plotSubnet
#'
#' This function allows you to get the adjacency matrix from a graph object.
#' @param pathway A graph object.
#' @param directed Set to TRUE if you want the graph to be directed.
#' @keywords adjacency matrix
#' @export
#' @examples
#' #get data
#' plotSubnet()

plotSubnet<-function(path,liste)
{
  subs<-as.character(liste)
  # subs<-translateGeneID2KEGGID(subs, organism="hsa")

  col<-c(rep("red",length(subs)))
  nAttrs<-list()
  nAttrs$fillcolor <- col
  nAttrs$height <- nAttrs$width <- rep("0.4", length(nodes(path)))
  names(nAttrs$width)<-names(nAttrs$height)<-nodes(path)
  names(nAttrs$fillcolor)<-c(subs)

  return(plot(path, y="neato", nodeAttrs = nAttrs))
}
