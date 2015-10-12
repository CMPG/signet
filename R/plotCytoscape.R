#' Plot the results in Cytoscape
#'
#' This function plots the output of the algorithm in Cytoscape.
#' CytoscapeRPC plugin must be installed and launched in Cytoscape (< 3.0)
#' before running this function.
#'
#' @param pathway A \verb{graphNEL} object.
#' @param outAlgo Simulated annealing output.
#' @param layout Layout in Cytoscape.
#'
#' @keywords plot, Cytoscape
#' @return A plot in Cytoscape.
#' @examples
#' require(signet)
#' data(keggPathways)

plotCytoscape <- function(pathway,outAlgo,layout="jgraph-spring")
{
  requireNamespace("RCytoscape",quietly=TRUE)

  # Nodes and edges attributes
  pathway <- initNodeAttribute(graph=pathway, attribute.name='subNet',
                               attribute.type='char',default.value='false')
  pathway <- initNodeAttribute(graph=pathway, attribute.name='z',
                               attribute.type='numeric',default.value=0)
  pathway <- initEdgeAttribute(graph=pathway, attribute.name='edgeType',
                               attribute.type='char',default.value='undefined')
  pathway <- initEdgeAttribute(graph=pathway, attribute.name='weight',
                               attribute.type='char',default.value='undefined')

  # add values corresponding to attributes (gene scores, ...)
  glist<-as.character(outAlgo$table[outAlgo$table$state,]$gene)
  for(gene in glist) nodeData(pathway,gene,"subNet")<-"true"
  glistz<-as.character(outAlgo$table$gene)
  for(gene in glistz) nodeData(pathway,gene,"z")<-outAlgo$table[outAlgo$table$gene==gene,]$score
  for(gene in glistz[!glistz%in%glist]) nodeData(pathway,gene,"subNet")<-"false"

  #create the network in Cytoscape
  cw <- new.CytoscapeWindow(p@title, graph=pathway,overwriteWindow=TRUE)
  displayGraph(cw)

  # setNodeShapeRule(cw,"subNet",c("true","false"),c('ellipse','ellipse'))

  #node colored in red if part of the subnetwork
  setNodeColorRule(cw,"subNet",c("true","false"),
                   c('#ff0000','#ffffff'),mode='lookup')

  #size of the node proportional to the gene score
  setNodeSizeRule(cw,"z",c(0,4),c(30,60),
                  mode= 'interpolate')

  setDefaultBackgroundColor(cw,'#ffffff') #white background
  setDefaultEdgeColor(cw,'#bdbdbd') #grey edge

  layoutNetwork(cw,layout) #update the network in Cytoscape
}
