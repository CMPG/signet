#' Undocumented functions
#' @export
#' @keywords internal

#returns summary statistics for all graphs in a database
graphSummary<-function(allgraphs){
  require(igraph)
  output<-matrix(NA,ncol=12,nrow=length(allgraphs))
  for(i in 1:length(allgraphs)){

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
  for(i in 2:12){
    output[,i]<-as.numeric(as.character(output[,i]))
  }
  return(output)
}


#' @export
#' @keywords internal
computeScore <- function(signetObject, score = "ideker")
{
  activeNet <- signetObject$network[signetObject$network$active, ]

  if (score == "mean") {
    subnetStat <- mean(activeNet$score)
  }
  if (score == "sum") {
    subnetStat <- sum(activeNet$score)
  }
  if (score == "ideker") {
    k <- length(activeNet$score)
    subnetStat <- (1/sqrt(k))*sum(activeNet$score)
  }
  if (score == "delta") {
    subnetStat <- mean(activeNet$score)-mean(tail(sort(activeNet$score),5))
  }
  return(subnetStat)
}



#' @export
#' @keywords internal
getAdjacencyMatrix<-function(pathway,
                             directed = FALSE,
                             selfLoops=FALSE)
{
  requireNamespace("graph",quietly=TRUE)
  x<-graph::edges(pathway)
  GList <- names(x)
  adjMatrix <- matrix(0,length(GList),length(GList),dimnames=list(GList,GList))

  if(!directed) {
    sapply(GList,function(y){
      adjMatrix[y,x[[y]]]<<-1
      adjMatrix[x[[y]],y]<<-1
    })
  } else {
    sapply(GList,function(y){
      adjMatrix[y,x[[y]]]<<-1
    })
  }

  if(!selfLoops) diag(adjMatrix)<-0

  rownames(adjMatrix)<-colnames(adjMatrix)<-GList

  return(adjMatrix)
}


#' @export
#' @keywords internal
#returns graphs
filterGraphs<-function(allgraphs,gsummary,nodes_min = 0,density_max = 1){
  filter<-which(gsummary$nodes >= nodes_min & gsummary$density <= density_max)
  return(allgraphs[filter])
}

#' @export
#' @keywords internal
#Output
returnTable<-function(outputSignet,
                      pDatabase,
                      nullDistribution,
                      density_cuts=0,
                      database = NA,
                      correction=TRUE)
{
  subnetSize<-unlist(lapply(outputSignet,function(x) {
    if(!is.null(x)){ stat<-x$subnet_size} else { stat<-NA}; return(stat)
  }))
  netSize<-unlist(lapply(outputSignet,function(x) {
    if(!is.null(x)){ stat<-length(x$network$gene)} else { stat<-NA}; return(stat)
  }))
  subnetScore<-unlist(lapply(outputSignet,function(x) {
    if(!is.null(x)){ stat<-x$subnet_score}else { stat<-NA}; return(stat)
  }))

  if(length(density_cuts) == 1){
    pvalues<-unlist(lapply(subnetScore,function(x) {
      stat<-mean(nullDistribution>x,na.rm=TRUE);if(is.null(stat)) stat<-NA; return(stat)
    }))
    pvalues[subnetSize<2]<-NA
  } else {
    X<-graphSummary(pDatabase)
    pvalues<-c()

    for(i in 1:(length(density_cuts)-1)) {
      ids<-which(X$density > density_cuts[i] & X$density < density_cuts[i+1] & X$nodes >= 10)

      for (ii in ids) {
        pvalues[ii]<-mean(nullDistribution[[i]]>subnetScore[ii],na.rm=TRUE);
        if(is.null(pvalues[ii])) pvalues[ii]<-NA;
      }
    }
  }

  pathwaysNames<-names(pDatabase)
  datab<-rep(database,length(pathwaysNames))
  if(correction){
    require(qvalue)
    qvalues<-qvalue(pvalues)$qvalues
  }
  geneListEntrez<-unlist(lapply(outputSignet,function(x) {
    stat<-paste(x$network[x$network$active,]$gene,collapse=";");
    if(is.null(stat)) stat<-NA; return(stat)
  }))

  if(correction){
    out<-data.frame(datab,pathwaysNames,netSize,subnetSize,subnetScore,
                    pvalues,qvalues,geneListEntrez)
  } else {
    out<-data.frame(datab,pathwaysNames,netSize,subnetSize,subnetScore,
                    pvalues,geneListEntrez)
  }

  return(out)
}

#' @export
#' @keywords internal
overlapMatrix<-function(signetObject){
  overlap<-matrix(0,nrow=length(signetObject),ncol=length(signetObject))
  for(i in 1:length(signetObject))
  {
    for(j in 1:length(signetObject))
    {
      gi<-signetObject[[i]]$network[signetObject[[i]]$network$active,]$gene
      gj<-signetObject[[j]]$network[signetObject[[j]]$network$active,]$gene
      inter<-length(intersect(gi,gj))
      un<-length(union(gi,gj))
      if(!is.null(inter)&!is.null(un)) overlap[i,j]<-inter/un
      else overlap[i,j]<-NA
    }
    print(i)
  }
  return(overlap)
}

#' @export
#' @keywords internal
correctOverlap<-function(overlapMatrix,
                         signetTable,
                         threshold,
                         pvalue = 0.01,
                         cex = 1,
                         graphics = TRUE) {

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
  # abline(h=c(100*(1-threshold)),col=2,lty=2,lwd=2)
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
      print(signetTable[which(signetTable$group==signetTable[i,]$group),1])
      minp<-min(signetTable[which(signetTable$group==signetTable[i,]$group),]$pvalues)
      print(minp)

      if(minp==signetTable[i,]$pvalues) {
        print(qvalmax$q[which(abs(qvalmax$pvalues-minp) == min(abs(qvalmax$pvalues-minp)))])
        signetTable[i,]$qvalues<-qvalmax$q[which(abs(qvalmax$pvalues-minp) == min(abs(qvalmax$pvalues-minp)))][1]
      }
    }
  }

  return(signetTable)
}

#' @export
#' @keywords internal
createSignetObject <- function(pathway, scores, iterations, minimumSize) {
  threshold <- minimumSize
  if (length(graph::nodes(pathway))==0) {
    X<-NA
  } else {
    TA<-unlist(lapply(graph::connComp(pathway),length))
    maxi<-which(TA==max(TA))
    CC<-unlist(graph::connComp(pathway)[maxi])
    X<-0
  }

  if (is.na(X)) {
    stop("\r  No connected component of size greater than 5.")
  }

  names(scores)[[1]]<-"gene";names(scores)[[2]]<-"score"
  scores<-scores[scores$gene %in% graph::nodes(pathway),]

  #to avoid duplicated nodes:
  X<-tapply(scores[!is.na(scores$score),]$score,scores[!is.na(scores$score),]$gene,max,na.rm=TRUE)
  X2<-data.frame(gene=as.numeric(names(X)),score=X)
  scores<-rbind(X2,scores[is.na(scores$score) & !scores$gene %in% X2$gene,])

  connected_comp<-graph::subGraph(as.character(scores[!is.na(scores$score),]$gene),pathway)
  connected_comp<-graph::subGraph(CC[CC%in%as.character(scores[!is.na(scores$score),]$gene)],connected_comp)
  connected_comp<-graph::ugraph(connected_comp)

  nnodes<-dim(scores)[1]
  subnet_score <- NA
  subnet_size <- NA
  subnet_genes <- NA
  network<-data.frame(
    gene=scores$gene,
    score=scores$score,
    active=rep(FALSE,nnodes)
  )

  if(max(TA) < minimumSize){
    object<-list(pathway=pathway,
                 connected_comp=connected_comp,
                 network=network,
                 subnet_score=subnet_score,
                 subnet_size=subnet_size,
                 subnet_genes=subnet_genes,
                 p.value=NA)
    class(object)<-"signet"
    return(object)
  } else {
    simulated_annealing<-data.frame(
      temperature=temperatureFunction(iterations),
      size_evolution=rep(NA,iterations),
      score_evolution=rep(NA,iterations)
    )

    object<-list(pathway=pathway,
                 connected_comp=connected_comp,
                 network=network,
                 subnet_score=subnet_score,
                 subnet_size=subnet_size,
                 subnet_genes=subnet_genes,
                 p.value=NA,
                 simulated_annealing=simulated_annealing)
    class(object)<-"signet"
    return(object)
  }
}

#' @export
#' @keywords internal
print.signet <- summary.signet <- function(object) {
  cat("High-scoring subnetwork found with simulated annealing\n")
  cat(paste("Subnetwork score: ",round(object$subnet_score,digits=4),"\n",sep=""))
  cat(paste("Subnetwork size: ",object$subnet_size,"\n",sep=""))
  cat(paste("Genes in subnetwork: ",paste(object$subnet_genes,collapse=" "),"\n",sep=""))
}

#' @export
#' @keywords internal
summary.signetList <- function(object) {
  signet_table <- data.frame(pathway=names(object),
                             net.size=unlist(lapply(object,function(x)dim(x$network)[1])),
                             subnet.size=unlist(lapply(object,function(x)x$subnet_size)),
                             subnet.score=unlist(lapply(object,function(x)x$subnet_score)),
                             p.value=unlist(lapply(object,function(x)x$p.value)),
                             subnet.genes=unlist(lapply(object,function(x)paste(x$subnet_genes,collapse=" "))))
  rownames(signet_table)<-NULL
  return(signet_table)
}

#' @export
#' @keywords internal
plot.signet <- function(object) {
  m <- rbind(c(0,1,1,1,1,0),c(2,2,2,3,3,3))
  layout(m)

  glist<-object$network[object$network$active,]$gene

  subs<-as.character(glist)
  subg<-graph::subGraph(subs,object$connected_comp)

  col<-c(rep("red",length(subs)))
  nAttrs<-list()
  nAttrs$fillcolor <- col
  nAttrs$height <- rep("0.6", length(graph::nodes(object$connected_comp)))
  nAttrs$width <- rep("0.6", length(graph::nodes(object$connected_comp)))
  nAttrs$color <- rep("darkgrey", length(graph::nodes(object$connected_comp)))

  names(nAttrs$color)<-names(nAttrs$width)<-names(nAttrs$height)<-graph::nodes(object$connected_comp)
  names(nAttrs$fillcolor)<-c(subs)

  eAttrs<-list()
  eAttrs$color <- rep("grey",length(graph::edgeNames(object$connected_comp)))
  names(eAttrs$color)<-graph::edgeNames(object$connected_comp)

  graph::plot(object$connected_comp, y="neato", nodeAttrs = nAttrs,
              edgeAttrs = eAttrs)

  x <- 1:length(object$simulated_annealing$temperature)
  y <- object$simulated_annealing$size_evolution
  z <- object$simulated_annealing$score_evolution
  par(mar = c(5, 5, 5, 5))  # Leave space for z axis

  plot(x, z, type = "l",lwd = 1,cex = 0.2,pch = 16,
       col = "firebrick",ylab = "Subnetwork score",
       xlab = "Iterations") # first plot

  par(new = TRUE)
  plot(x, object$simulated_annealing$temperature, type = "l", lwd=1,
       xlab = "", ylab = "", axes=F, lty=2,
       col = "grey")
  axis(side = 4, at = pretty(range(object$simulated_annealing$temperature)))
  mtext("Temperature", side = 4, line = 3, cex=0.7)

  plot(x, y, type = "l", cex = 0.5,
       xlab = "Iterations", ylab = "Subnetwork size",
       pch=16,lwd=1,col="dodgerblue")

  par(new=TRUE)
  plot(x,object$simulated_annealing$temperature, type="l",lwd=1,
       xlab="", ylab="",axes=F,lty=2,
       col="grey")
  axis(side=4, at = pretty(range(object$simulated_annealing$temperature)))
  mtext("Temperature", side=4, line=3, cex=0.7)

  par(mfrow = c(1,1))

}

#' @export
#' @keywords internal
adjacencyMatrixToList <- function(adjMatrix) {
  adjList <- apply(adjMatrix, 1, function(x) return(names(x[x==1])))
  return(adjList)
}

#' @export
#' @keywords internal
plot.Cytoscape<-function(object,graphlist,results,threshold){
  clust<-object[which(object$pvalue<threshold),]
  #get number of occurences of a gene (color)
  glist<-strsplit(as.character(clust[1,]$geneListEntrez), ";")[[1]]
  ID<-which(names(graphlist)==clust$pathwaysNames[1])
  GRAPH<-subGraph(glist,graphlist[ID][[1]])
  GRAPH<-graph_from_graphnel(GRAPH)
  for(i in 2:dim(clust)[1]){
    glist2<-strsplit(as.character(clust[i,]$geneListEntrez), ";")[[1]]
    glistInGraph<-glist2[glist2%in%nodes(graphlist[[as.character(clust$pathwaysNames[i])]])]
    subG2<-subGraph(glistInGraph,graphlist[[as.character(clust$pathwaysNames[i])]])
    subG2<-graph_from_graphnel(subG2)

    GRAPH<-graph.union(GRAPH,subG2)
    E(GRAPH)$weight_1[is.na(E(GRAPH)$weight_1)]<-0
    E(GRAPH)$weight_2[is.na(E(GRAPH)$weight_2)]<-0

    E(GRAPH)$weight <- E(GRAPH)$weight_1 + E(GRAPH)$weight_2
  }
  Glist<-names(GRAPH[1])

  #create TABLE with gene and score and occurences
  TABLE<-data.frame(gene=Glist,score=rep(NA,length(Glist)),occurence=rep(0,length(Glist)))
  for(i in 1:dim(clust)[1]){ #scores
    ID<-which(names(graphlist)==clust$pathwaysNames[i])
    Gscor<-results[ID][[1]]$network[results[ID][[1]]$network$active,]$score
    Gnam<-results[ID][[1]]$network[results[ID][[1]]$network$active,]$gene
    for(ii in 1:length(Gscor)){
      if(sum(TABLE$gene==as.character(Gnam[ii]))>0){
        TABLE[TABLE$gene==as.character(Gnam[ii]),]$score<-Gscor[ii]
      }
    }
  }

  #nodes occurences
  for(i in 1:dim(clust)[1]){ #scores
    ID<-which(names(graphlist)==clust$pathwaysNames[i])
    Gnam<-results[ID][[1]]$network[results[ID][[1]]$network$active,]$gene
    for(ii in 1:length(Gnam)){
      if(sum(TABLE$gene==as.character(Gnam[ii]))>0){
        TABLE[TABLE$gene==as.character(Gnam[ii]),]$occurence<-TABLE[TABLE$gene==as.character(Gnam[ii]),]$occurence+1
      }
    }
  }

  #send to cytoscape
  library(RCytoscape);library(biomaRt)
  graph<-ugraph(as_graphnel(GRAPH))
  ex<-nodes(graph)

  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene2genomeEx <- getBM(values = ex, filters = "entrezgene", mart = mart, attributes = c("hgnc_symbol","entrezgene"))

  for(i in 1:length(nodes(graph))){
    # nodes(graph)[i]<-gene2genomeEx[which(gene2genomeEx$entrezgene==nodes(graph)[i]),]$hgnc_symbol
    interlist<-gene2genomeEx[which(gene2genomeEx$entrezgene==TABLE$gene[i]),]$hgnc_symbol
    if(length(interlist)[1]>1) {
      nodes(graph)[i]<-interlist[nchar(interlist)>0][1]
    } else {
      nodes(graph)[i]<-interlist
    }
  }
  TABLE$gene<-as.character(TABLE$gene)
  for(i in 1:length(TABLE$gene)){
    interlist<-gene2genomeEx[which(gene2genomeEx$entrezgene==TABLE$gene[i]),]$hgnc_symbol
    if(length(interlist)[1]>1) {
      TABLE$gene[i]<-interlist[nchar(interlist)>0][1]
    } else {
      TABLE$gene[i]<-interlist
    }
  }

  graph <- initNodeAttribute(graph=graph,  attribute.name='moleculeType',attribute.type='char',default.value='false')
  graph <- initNodeAttribute(graph=graph,  attribute.name='occurence',attribute.type='numeric',default.value=1)
  graph <- initNodeAttribute(graph=graph,  attribute.name='score',attribute.type='numeric',default.value=0)

  graph <- initEdgeAttribute(graph=graph, attribute.name='edgeType',attribute.type='char',default.value='undefined')
  graph <- initEdgeAttribute(graph=graph, attribute.name='weight',attribute.type='char',default.value=0)

  glistz<-nodes(graph)
  for(gene in glistz) nodeData(graph,gene,"score")<-TABLE[TABLE$gene==gene,]$score
  for(gene in glistz) nodeData(graph,gene,"occurence")<-TABLE[TABLE$gene==gene,]$occurence


  cw <- new.CytoscapeWindow("signet", graph=graph,overwriteWindow=TRUE)

  displayGraph(cw)

  setNodeColorRule(cw,"score",c(0,1),c('#ffffff','#b22222'),mode= 'interpolate')
  setNodeSizeRule(cw,"score",c(0,1),c(40,70),mode= 'interpolate')
  # setNodeBorderWidthRule(cw,"occurence",c(1:max(TABLE$occurence)),c(1,max(TABLE$occurence)),1)
  setEdgeLineWidthRule(cw,"weight",c(1,10),c(1,50))
  setDefaultBackgroundColor(cw,'#ffffff')
  setDefaultEdgeColor(cw,'#bdbdbd')
  layoutNetwork(cw,"jgraph-spring")

}

#' @export
#' @keywords internal
plot.single.Cytoscape<-function(object,pnumber,graphlist,results){
  clust<-object[pnumber,]
  #get number of occurences of a gene (color)
  glist<-strsplit(as.character(clust[1,]$geneListEntrez), ";")[[1]]
  ID<-pnumber
  GRAPH<-graphlist[ID][[1]]#subGraph(glist,graphlist[ID][[1]])
  GRAPH<-graph_from_graphnel(GRAPH)
  Glist<-names(GRAPH[1])

  #create TABLE with gene and score and occurences
  TABLE<-data.frame(gene=Glist,score=rep(NA,length(Glist)),occurence=rep(0,length(Glist)))
  # for(i in 1:dim(clust)[1]){ #scores
  ID<-which(names(graphlist)==clust$pathwaysNames)
  Gscor<-results[ID][[1]]$network$score
  Gnam<-results[ID][[1]]$network$gene
  for(ii in 1:length(Gscor)){
    if(sum(TABLE$gene==as.character(Gnam[ii]))>0){
      TABLE[TABLE$gene==as.character(Gnam[ii]),]$score<-Gscor[ii]
    }
  }
  # }
  #nodes occurences
  Gnam<-results[ID][[1]]$network[results[ID][[1]]$network$active,]$gene
  TABLE[TABLE$gene %in% as.character(Gnam),]$occurence<-1

  #send to cytoscape
  library(RCytoscape);library(biomaRt)
  graph<-ugraph(as_graphnel(GRAPH))
  ex<-nodes(graph)

  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene2genomeEx <- getBM(values = ex, filters = "entrezgene", mart = mart, attributes = c("hgnc_symbol","entrezgene"))

  for(i in 1:length(nodes(graph))){
    nodes(graph)[i]<-gene2genomeEx[which(gene2genomeEx$entrezgene==nodes(graph)[i]),]$hgnc_symbol
  }
  TABLE$gene<-as.character(TABLE$gene)
  for(i in 1:length(TABLE$gene)){
    TABLE$gene[i]<-gene2genomeEx[which(gene2genomeEx$entrezgene==TABLE$gene[i]),]$hgnc_symbol
  }

  TABLE$score[is.na(TABLE$score)]<-0
  graph <- initNodeAttribute(graph=graph,  attribute.name='moleculeType',attribute.type='char',default.value='false')
  graph <- initNodeAttribute(graph=graph,  attribute.name='occurence',attribute.type='numeric',default.value=1)
  graph <- initNodeAttribute(graph=graph,  attribute.name='score',attribute.type='numeric',default.value=0)

  graph <- initEdgeAttribute(graph=graph, attribute.name='edgeType',attribute.type='char',default.value='undefined')
  graph <- initEdgeAttribute(graph=graph, attribute.name='weight',attribute.type='char',default.value=0)

  glistz<-nodes(graph)
  for(gene in glistz) nodeData(graph,gene,"score")<-TABLE[TABLE$gene==gene,]$score
  for(gene in glistz) nodeData(graph,gene,"occurence")<-TABLE[TABLE$gene==gene,]$occurence


  cw <- new.CytoscapeWindow("Coverage", graph=graph,overwriteWindow=TRUE)

  displayGraph(cw)

  setNodeColorRule(cw,"occurence",c(0,1),c('#ffffff','#b22222'),mode= 'interpolate')
  setNodeSizeRule(cw,"score",c(min(TABLE$score),max(TABLE$score)),c(40,70),mode= 'interpolate')
  # setNodeBorderWidthRule(cw,"occurence",c(1:max(TABLE$occurence)),c(1,max(TABLE$occurence)),1)
  setEdgeLineWidthRule(cw,"weight",c(1,10),c(1,50))
  setDefaultBackgroundColor(cw,'#ffffff')
  setDefaultEdgeColor(cw,'#bdbdbd')
  layoutNetwork(cw,"jgraph-spring")

}

#' @export
#' @keywords internal
lapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}
