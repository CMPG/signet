#' @keywords internal
overlapMatrix<-function(signetObject) {

  overlap<-matrix(0,nrow=length(signetObject),ncol=length(signetObject))

  for(i in 1:length(signetObject)) {

    for(j in 1:length(signetObject)) {

      gi<-signetObject[[i]]$network[signetObject[[i]]$network$active,]$gene
      gj<-signetObject[[j]]$network[signetObject[[j]]$network$active,]$gene
      inter<-length(intersect(gi,gj))
      un<-length(union(gi,gj))

      if(!is.null(inter)&!is.null(un)) overlap[i,j]<-inter/un
      else overlap[i,j]<-NA

    }

  }

  return(overlap)

}



#' @keywords internal
correctOverlap<-function(overlapMatrix,
                         signetTable,
                         threshold,
                         pvalue = 0.01,
                         cex = 1,
                         graphics = TRUE) {

  requireNamespace("stats", quietly = TRUE)

  if(graphics) {
    torm <- (!is.na(signetTable$pvalues) & signetTable$pvalues < pvalue)
    overlap2 <- overlapMatrix[c(torm),c(torm)]
    d <- as.dist(100*(1-overlap2))
    h <- hclust(d)

    plot(h, main = "Overlapping correction", sub = "", xlab="",
         ylab="Overlapping percentage between subnetworks",
         axes = FALSE, hang = -1,
         labels=paste(signetTable$pathwaysNames[torm]," (p = ",
                      round(signetTable$pvalues[torm],digits=3),")",sep=""),
         cex=cex)
    lines(x = c(0,0), y = c(0,100), type = "n") # force extension of y axis
    axis(side = 2, at = seq(0,100,10), labels = seq(100,0,-10))
  }

  torm <- !is.na(signetTable$pvalues)
  overlap3 <- overlapMatrix[c(torm),c(torm)]
  d2 <- as.dist(100*(1-overlap3))
  h2 <- hclust(d2)

  clust <- as.factor(cutree(h2,h=100*(1-threshold)))
  group<-torm
  group[group]<-clust
  group[!group]<-NA
  signetTable$group<-group

  pvalmax<-tapply(signetTable$pvalues,signetTable$group,max)

  require(qvalue)

  qvalmax<-qvalue(pvalmax)
  signetTable$qvalues<-rep(NA,length(signetTable$pvalues))

  for(i in 1:length(signetTable$pvalues)) {

    if(!is.na(signetTable[i,]$pvalues)) {

      minp<-min(signetTable[which(signetTable$group ==
                                    signetTable[i,]$group),]$pvalues)

      if(minp==signetTable[i,]$pvalues) {

        print(qvalmax$q[which(abs(qvalmax$pvalues-minp) ==
                                min(abs(qvalmax$pvalues-minp)))])

        signetTable[i,]$qvalues <-
          qvalmax$q[which(abs(qvalmax$pvalues-minp) ==
                            min(abs(qvalmax$pvalues-minp)))][1]
      }
    }
  }

  return(signetTable)
}

#' @keywords internal

#returns summary statistics for all graphs in a database
graphSummary<-function(allgraphs) {

  require(igraph)

  output<-matrix(NA,ncol=12,nrow=length(allgraphs))

  for(i in 1:length(allgraphs)) {

    g<-igraph.from.graphNEL(allgraphs[[i]], name = TRUE, weight = TRUE,
                            unlist.attrs = TRUE)

    comp<-components(g)
    g2<- induced.subgraph(g,names(comp$membership[comp$membership==1]))

    bla<-c(length(V(g)),
           length(E(g)),
           graph.density(g),
           diameter(g),
           clusters(g)$no,
           length(V(g2)),
           length(E(g2)),
           graph.density(g2),
           diameter(g2),
           edge.connectivity(g2),
           graph.adhesion(g2))

    output[i,1]<-names(allgraphs[i])
    output[i,2:12]<-as.numeric(as.character(bla))

  }

  output<-as.data.frame(output)

  colnames(output)<-c("name",
                      "nodes",
                      "edges",
                      "density",
                      "diameter",
                      "islands",
                      "LCC_nodes",
                      "LCC_edges",
                      "LCC_density",
                      "LCC_diameter",
                      "LCC_edges_connectivity",
                      "LCC_adhesion")

  for(i in 2:12) {

    output[,i]<-as.numeric(as.character(output[,i]))

  }

  return(output)

}


#' @keywords internal
plot.Cytoscape<-function(object,graphlist,results,threshold) {

  require(igraph)
  object<-object[!is.nan(object$p.value),]
  object<-object[!is.na(object$p.value),]
  clust<-object[which(object$p.value<threshold),]

  #get number of occurences of a gene (color)
  glist<-strsplit(as.character(clust[1,]$subnet.genes), " ")[[1]]
  ID<-which(names(graphlist)==clust$pathway[1])
  GRAPH<-graph::subGraph(glist,graphlist[ID][[1]])
  GRAPH<-igraph::graph_from_graphnel(GRAPH)

  for(i in 2:dim(clust)[1]) {

    glist2<-strsplit(as.character(clust[i,]$subnet.genes), " ")[[1]]
    glistInGraph<-glist2[
      glist2 %in% graph::nodes(graphlist[[as.character(clust$pathway[i])]])]

    subG2<-subGraph(glistInGraph,graphlist[[as.character(clust$pathway[i])]])
    subG2<-igraph::graph_from_graphnel(subG2)

    GRAPH<-igraph::graph.union(GRAPH,subG2)
    E(GRAPH)$weight_1[is.na(E(GRAPH)$weight_1)]<-0
    E(GRAPH)$weight_2[is.na(E(GRAPH)$weight_2)]<-0

    E(GRAPH)$weight <- E(GRAPH)$weight_1 + E(GRAPH)$weight_2
  }
  Glist<-names(GRAPH[1])

  #create TABLE with gene and score and occurences
  TABLE<-data.frame(gene=Glist,
                    score=rep(NA,length(Glist)),
                    occurence=rep(0,length(Glist)))

  for(i in 1:dim(clust)[1]) { #scores

    ID<-which(names(graphlist)==clust$pathway[i])
    Gscor<-results[ID][[1]]$network[results[ID][[1]]$network$active,]$score
    Gnam<-results[ID][[1]]$network[results[ID][[1]]$network$active,]$gene

    for(ii in 1:length(Gscor)) {
      if(sum(TABLE$gene==as.character(Gnam[ii]))>0){
        TABLE[TABLE$gene==as.character(Gnam[ii]),]$score<-Gscor[ii]
      }
    }

  }

  #nodes occurences
  for(i in 1:dim(clust)[1]){ #scores

    ID<-which(names(graphlist)==clust$pathway[i])
    Gnam<-results[ID][[1]]$network[results[ID][[1]]$network$active,]$gene

    for(ii in 1:length(Gnam)) {

      if(sum(TABLE$gene==as.character(Gnam[ii]))>0) {
        TABLE[TABLE$gene==as.character(Gnam[ii]),]$occurence <-
          TABLE[TABLE$gene==as.character(Gnam[ii]),]$occurence+1
      }

    }

  }

  #send to cytoscape
  library(RCytoscape);library(biomaRt)
  graph<-ugraph(as_graphnel(GRAPH))
  ex<-graph::nodes(graph)

  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene2genomeEx <- getBM(values = ex,
                         filters = "entrezgene",
                         mart = mart,
                         attributes = c("hgnc_symbol","entrezgene"))

  for(i in 1:length(graph::nodes(graph))) {

    interlist<-gene2genomeEx[which(gene2genomeEx$entrezgene ==
                                     TABLE$gene[i]),]$hgnc_symbol

    if(length(interlist)[1]>1) {
      graph::nodes(graph)[i]<-interlist[nchar(interlist)>0][1]
    } else {
      graph::nodes(graph)[i]<-interlist
    }
  }

  TABLE$gene<-as.character(TABLE$gene)

  for(i in 1:length(TABLE$gene)) {

    interlist<-gene2genomeEx[which(gene2genomeEx$entrezgene ==
                                     TABLE$gene[i]),]$hgnc_symbol

    if(length(interlist)[1]>1) {
      TABLE$gene[i]<-interlist[nchar(interlist)>0][1]
    } else {
      TABLE$gene[i]<-interlist
    }

  }

  graph <- initNodeAttribute(graph=graph,
                             attribute.name='moleculeType',
                             attribute.type='char',
                             default.value='false')
  graph <- initNodeAttribute(graph=graph,
                             attribute.name='occurence',
                             attribute.type='numeric',
                             default.value=1)
  graph <- initNodeAttribute(graph=graph,
                             attribute.name='score',
                             attribute.type='numeric',
                             default.value=0)

  graph <- initEdgeAttribute(graph=graph,
                             attribute.name='edgeType',
                             attribute.type='char',
                             default.value='undefined')
  graph <- initEdgeAttribute(graph=graph,
                             attribute.name='weight',
                             attribute.type='char',
                             default.value=0)

  glistz<-graph::nodes(graph)
  for(gene in glistz) nodeData(graph,gene,"score")<-TABLE[TABLE$gene==gene,]$score
  for(gene in glistz) nodeData(graph,gene,"occurence")<-TABLE[TABLE$gene==gene,]$occurence


  cw <- new.CytoscapeWindow("signet", graph=graph,overwriteWindow=TRUE)

  displayGraph(cw)

  setNodeColorRule(cw,"score",c(0,1),c('#ffffff','#b22222'),mode= 'interpolate')
  setNodeSizeRule(cw,"score",c(0,1),c(40,70),mode= 'interpolate')

  setEdgeLineWidthRule(cw,"weight",c(1,10),c(1,50))
  setDefaultBackgroundColor(cw,'#ffffff')
  setDefaultEdgeColor(cw,'#bdbdbd')
  layoutNetwork(cw,"jgraph-spring")

}


#' @keywords internal
plot.single.Cytoscape<-function(object,pnumber,graphlist,results) {

  clust<-object[pnumber,]

  #get number of occurences of a gene (color)
  glist<-strsplit(as.character(clust[1,]$geneListEntrez), ";")[[1]]
  ID<-pnumber
  GRAPH<-graphlist[ID][[1]]
  GRAPH<-graph_from_graphnel(GRAPH)
  Glist<-names(GRAPH[1])

  #create TABLE with gene and score and occurences
  TABLE<-data.frame(gene=Glist,
                    score=rep(NA,length(Glist)),occurence=rep(0,length(Glist)))

  ID<-which(names(graphlist)==clust$pathwaysNames)
  Gscor<-results[ID][[1]]$network$score
  Gnam<-results[ID][[1]]$network$gene

  for(ii in 1:length(Gscor)) {
    if(sum(TABLE$gene==as.character(Gnam[ii]))>0){
      TABLE[TABLE$gene==as.character(Gnam[ii]),]$score<-Gscor[ii]
    }
  }

  #nodes occurences
  Gnam<-results[ID][[1]]$network[results[ID][[1]]$network$active,]$gene
  TABLE[TABLE$gene %in% as.character(Gnam),]$occurence<-1

  #send to cytoscape
  library(RCytoscape);library(biomaRt)
  graph<-ugraph(as_graphnel(GRAPH))
  ex<-nodes(graph)

  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene2genomeEx <- getBM(values = ex,
                         filters = "entrezgene",
                         mart = mart,
                         attributes = c("hgnc_symbol","entrezgene"))

  for(i in 1:length(nodes(graph))) {
    nodes(graph)[i]<-gene2genomeEx[which(gene2genomeEx$entrezgene ==
                                           nodes(graph)[i]),]$hgnc_symbol
  }

  TABLE$gene<-as.character(TABLE$gene)

  for(i in 1:length(TABLE$gene)) {
    TABLE$gene[i]<-gene2genomeEx[which(gene2genomeEx$entrezgene ==
                                         TABLE$gene[i]),]$hgnc_symbol
  }

  TABLE$score[is.na(TABLE$score)]<-0

  graph <- initNodeAttribute(graph=graph,
                             attribute.name='moleculeType',
                             attribute.type='char',
                             default.value='false')
  graph <- initNodeAttribute(graph=graph,
                             attribute.name='occurence',
                             attribute.type='numeric',
                             default.value=1)
  graph <- initNodeAttribute(graph=graph,
                             attribute.name='score',
                             attribute.type='numeric',
                             default.value=0)

  graph <- initEdgeAttribute(graph=graph,
                             attribute.name='edgeType',
                             attribute.type='char',
                             default.value='undefined')
  graph <- initEdgeAttribute(graph=graph,
                             attribute.name='weight',
                             attribute.type='char',
                             default.value=0)

  glistz<-nodes(graph)
  for(gene in glistz) {
    nodeData(graph,gene,"score") <- TABLE[TABLE$gene==gene,]$score
  }
  for(gene in glistz) {
    nodeData(graph,gene,"occurence") <- TABLE[TABLE$gene==gene,]$occurence
  }

  cw <- new.CytoscapeWindow("Coverage", graph=graph,overwriteWindow=TRUE)

  displayGraph(cw)

  setNodeColorRule(cw,"occurence",c(0,1),
                   c('#ffffff','#b22222'),mode= 'interpolate')
  setNodeSizeRule(cw,"score",c(min(TABLE$score),
                               max(TABLE$score)),c(40,70),mode= 'interpolate')

  setEdgeLineWidthRule(cw,"weight",c(1,10),c(1,50))
  setDefaultBackgroundColor(cw,'#ffffff')
  setDefaultEdgeColor(cw,'#bdbdbd')
  layoutNetwork(cw,"jgraph-spring")

}
