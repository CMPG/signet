#' Undocumented functions
#' @export

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

#returns graphs
filterGraphs<-function(allgraphs,gsummary,nodes_min = 0,density_max = 1){
  filter<-which(gsummary$nodes >= nodes_min & gsummary$density <= density_max)
  return(allgraphs[filter])
}


#O
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
    if(!is.null(x)){ stat<-x$subnet_score*(sqrt(x$subnet_size)/x$subnet_size)}else { stat<-NA}; return(stat)
  }))
  print(subnetScore)

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


overlapMatrix<-function(signetObject){
  overlap<-matrix(0,nrow=length(signetObject),ncol=length(signetObject))
  for(i in 1:length(signetObject))
  {
    for(j in 1:length(signetObject))
    {
      gi<-signetObject[[i]]$table[signetObject[[i]]$table$state,]$gene
      gj<-signetObject[[j]]$table[signetObject[[j]]$table$state,]$gene
      inter<-length(intersect(gi,gj))
      un<-length(union(gi,gj))
      if(!is.null(inter)&!is.null(un)) overlap[i,j]<-inter/un
      else overlap[i,j]<-NA
    }
    print(i)
  }
  return(overlap)
}

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
  abline(h=c(100*(1-threshold)),col=2,lty=2,lwd=2)
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



  for(i in 1:length(signetTable$pvalues)){


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

  #need ids of remaining subnetworks

  #compute qvalue using FDR method

  #add qvalues to table

  #plot histogram
  return(signetTable)
}

createSignetObject <- function(pathway, scores, iterations, minimumSize) {

  threshold <- minimumSize
  if (length(graph::nodes(pathway))==0) X<-NULL
  else  X<-unlist(graph::connComp(pathway)[!lapply(graph::connComp(pathway),
                                                         length)<threshold])
  if (is.null(X)) {
    stop("\r  No connected component of size greater than the specified threshold")
  } else {
    connected_comp<-graph::subGraph(X,pathway)
    connected_comp<-graph::ugraph(connected_comp)
  }

  names(scores)[[1]]<-"gene";names(scores)[[2]]<-"score"
  scores<-scores[scores$gene %in% graph::nodes(pathway),]

  connected_comp<-graph::subGraph(as.character(scores[!is.na(scores$score),]$gene),pathway)

  nnodes<-dim(scores)[1]
  subnet_score <- NA
  subnet_size <- NA
  subnet_genes <- rep(NA,nnodes)
  network<-data.frame(
    gene=scores$gene,
    score=scores$score,
    active=rep(FALSE,nnodes)
  )
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
               simulated_annealing=simulated_annealing)
  class(object)<-"signet"
  return(object)
}

print.signet <- summary.signet <- function(object) {
  cat("High-scoring subnetwork found with simulated annealing\n\n")
  cat(paste("  Subnetwork score: ",object$subnet_score,"\n",sep=""))
  cat(paste("  Subnetwork size: ",object$subnet_size,"\n",sep=""))
  cat(paste("  Genes in subnetwork: ",object$subnet_genes,"\n",sep=""))
}

plot.signet <- function(object) {
  m <- rbind(c(0,1,1,0),c(2,2,3,3))
  layout(m)

  glist<-object$network[object$network$active,]$genes

  subs<-as.character(glist)
  subg<-graph::subGraph(subs,object$connected_comp)

  col<-c(rep("red",length(subs)))
  nAttrs<-list()
  nAttrs$fillcolor <- col
  nAttrs$height <- nAttrs$width <- rep("0.4", length(graph::nodes(object$connected_comp)))
  names(nAttrs$width)<-names(nAttrs$height)<-graph::nodes(object$connected_comp)
  names(nAttrs$fillcolor)<-c(subs)

  graph::plot(object$connected_comp, y="neato", nodeAttrs = nAttrs,graph=list(overlap="scale"))

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
  mtext("Temperature", side = 4, line = 3)

  plot(x, y, type = "l", cex = 0.5,
       xlab = "Iterations", ylab = "Subnetwork size",
       pch=16,lwd=1,col="dodgerblue")

  par(new=TRUE)
  plot(x,object$simulated_annealing$temperature, type="l",lwd=1,
       xlab="", ylab="",axes=F,lty=2,
       col="grey")
  axis(side=4, at = pretty(range(object$simulated_annealing$temperature)))
  mtext("Temperature", side=4, line=3)

  par(mfrow = c(1,1))

}

adjacencyMatrixToList <- function(adjMatrix) {
  adjList <- apply(adjMatrix, 1, function(x) return(names(x[x==1])))
  return(adjList)
}

bfs <- function(graph, seeds) {

  # graph<-remove_loops(graph)
  V <- names(graph)
  N <- length(V)
  reachable <- seeds
  marks <- structure(rep(FALSE, N), names = V)

  while (length(seeds)) {

    s <- seeds[1]
    seeds <- seeds[-1]

    for (n in graph[[s]]) {
      if (!marks[[n]]) {
        seeds <- c(seeds, n)
        reachable <- c(reachable, n)
        marks[[n]] <- TRUE
      }
    }
  }

  reachable
}

remove_loops <- function(graph) {
  structure(
    lapply(names(graph), function(n) setdiff(graph[[n]], n)),
    names = names(graph)
  )
}

isArticulationPoint <- function(node, adjMatrix) {

  adjNodes <- names(which(adjMatrix[node, ] == 1))

  if(length(adjNodes) == 1) {
    isTrue <- FALSE
  } else {
    cond <- !(colnames(adjMatrix) %in% node)
    newMatrix <- adjMatrix[cond, cond]
    newList <- adjacencyMatrixToList(newMatrix)

    reachable <- list()
    for(i in adjNodes) {
      reachable[[i]] <- bfs(newList,adjNodes[adjNodes==i])
    }

    isTrue <- unlist(lapply(adjNodes,
                            function(x) {
                              sum(reachable[[x]] %in% adjNodes[adjNodes != x])
                              }))
    if(any(isTrue == 0)) {
      isTrue <- TRUE
    } else {
      isTrue <- FALSE
    }
  }
  return(isTrue)
}
