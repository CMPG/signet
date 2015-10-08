#' Search for a high scoring subnetwork
#'
#' A simulated annealing algorithm is used to research a high scoring subnetwork.
#' Implements also a test id a null distribution is provided.
#' @param pathway A graph object.
#' @param scores a data frame with gene list and associated scores
#' @param nullDist A data frame with three columns : k, mu, sigma. Can be obtained thanks to nulldistribution() function
#' @param iterations
#' @param temperature The temperature function parameter.
#' @param kmin The minimal size of a subnetwork.
#' @param directed If TRUE, considers the edges direction, i.e. cannot go back.
#' @param verbose If TRUE, displays text in the R console.
#' @keywords subnetwork, simulated annealing, heuristics, search algorithm
#' @export
#' @examples
#' data(keggPathways)
#' data(zScores)
#'
#' searchSubnet(keggPathways[[1]],zScores)

searchSubnet<-function(pathway,
                       scores,
                       nullDist,
                       iterations = 2000,
                       temperature = 0.995,
                       kmin = 3,
                       directed = FALSE,
                       verbose = TRUE,
                       animPlot = 0)
{
  #check for the package graph and load it
  if(sum(installed.packages()[,1]=="graph")==0)
  {
    stop("Package graph is not installed !")
  }
  else
  {
    suppressMessages(require(graph))
  }
  #check for arguments
  if(missing(pathway)) stop("Pathway is missing !")
  if(missing(scores)) stop("Gene scores are missing !")

  #remove conn components of size < kmin
  pathway<-subGraph(unlist(connComp(pathway)[!lapply(connComp(pathway),length)<kmin]),pathway)

  if(directed==FALSE) pathway<-ugraph(pathway)

  names(scores)[[1]]<-"gene";names(scores)[[2]]<-"score"
  scores<-scores[scores$gene%in%nodes(pathway),]

  if(missing(nullDist))
  {
    if(length(scores$gene)<100) max<-length(scores$gene)
    else max <- 100
    nullDist<-nullDistribution(scores = scores, kmin = 2, kmax = max,iterations=2000)
  }

  ### Initialization of the graph

  newpath<-subGraph(as.character(scores[!is.na(scores$score),]$gene),pathway) #remove missing values
  adjMatrix<-getAdjacencyMatrix(pathway=newpath) #get the adjacency matrix of the graph

  #initialize the gene scores list (gene names, scores, and state)
  workingTable <- scores[!is.na(scores$score),]
  workingTable$state <- rep(FALSE,length(workingTable$score))

  #initialize the temperature function
  temp<-temperatureFunction(iterations = iterations, param = temperature, burnin = 50)

  #c. select the initial subset of size kmin
  boundaries<-NULL
  geneSampled<-NULL
  for(ii in 1:kmin)
  {
    bu<-0
    gene<-NULL
    while(sum(bu)==0)
    {
      if(ii == 1){ gene<-sample(colnames(adjMatrix),1)}## 1. 1 gene au hasard
      else{ gene<-sample(boundaries,1)}
      bu<-adjMatrix[gene,]
    }
    boundaries<-c(boundaries,names(bu[which(bu!=0)]))

    boundaries<-boundaries[!boundaries %in% geneSampled]
    if(ii==1)
    {
      geneSampled<-gene
    }
    else
    {
      if(length(boundaries)>0) geneSampled<-c(geneSampled,gene)
    }
  }

  #toggle state in final list
  workingTable[which(workingTable$gene%in%geneSampled),]$state <- TRUE
  workingTable[which(!workingTable$gene%in%geneSampled),]$state <- FALSE
  sumStat<-mean(workingTable[workingTable$state,]$score)
  s<-(sumStat-nullDist[nullDist$k==kmin,]$mu)/nullDist[nullDist$k==kmin,]$sigma

  ###2. OPTIMISATION
  ######

  for(i in 1:iterations)
  {

    if(verbose == TRUE & (100*i/iterations)%%10==0)#displays the progession every 10%
    {
      cat('\r  Searching for a high-scoring subnetwork...',paste(10*(100*i/iterations)%/%10,"%",sep=""))
      flush.console()
    }

    ###
    activeNet<-as.character(workingTable[workingTable$state,]$gene)
    # subG<-subGraph(activeNet,newpath)
    # boundaries<-boundary(activeNet,newpath)
    bu<-adjMatrix[activeNet,]

    if(length(dim(bu))>0)
    {
      boundaries<-bu[apply(bu,1,sum)>0,]
      boundaries<-colnames(boundaries[,apply(boundaries,2,sum)>0])
      #remove boundaries already active
      boundaries<-boundaries[!boundaries %in% workingTable[workingTable$state,]$gene]
    }
    else boundaries<-NULL

    bla<-0
    while(bla!=1)
    {
      ge<-sample(activeNet,1)
      test<-subGraph(as.character(activeNet[activeNet!=ge]),newpath)

      bla<-length(connComp(test))
      if(bla==1) neighbours<-c(ge)
    }

    pro<-1
    if(length(activeNet)>2)
    {
      pro<-length(boundaries)/(length(boundaries)+length(activeNet))
    }

    sampBound<-FALSE
    if(pro>runif(1))
    {
      sampBound<-TRUE
    }

    if(length(activeNet)>=kmin & length(unique(c(neighbours,boundaries)))>0 & length(activeNet)<=100)
    {
      #Sample a new gene to toggle its state
      if(sampBound) newG<-as.character(sample(unique(c(neighbours,boundaries)),1))
      else newG<-neighbours
      workingTable[which(workingTable$gene==newG),]$state <- !workingTable[which(workingTable$gene==newG),]$state

      #Compute the subnet score
      sumStat<-mean(workingTable[workingTable$state,]$score)
      s2<-(sumStat-nullDist[nullDist$k==length(activeNet),]$mu)/nullDist[nullDist$k==length(activeNet),]$sigma

      #Keep or not the toggled gene

      if(s2<s)#If the score is weaker than the last one, we keep the gene toggled with a probability "prob"
      {
        prob<-exp((s2-s)/temp[i])
        if(prob>runif(1))
        {
          s<-s2
        }
        else
        {
          workingTable[which(workingTable$gene==newG),]$state<-!workingTable[which(workingTable$gene==newG),]$state
        }
      }
      else if(s2>s)#If the score is greater than the last one, we keep the gene toggled
      {
        s<-s2
      }
    }

    if(animPlot>0 & i<animPlot)
    {
      if(i==1)size<-NULL
      size<-c(size,length(workingTable[workingTable$state,]$gene))
      if(i==1)score<-NULL
      score<-c(score,s)

      par(mfrow=c(2,2))
      plotSubnet(newpath,workingTable[workingTable$state,]$gene)
      par(mar=c(5,5,5,5))
      plot(1:i,size,ylim=c(0,50),xlim=c(0,animPlot),
           ylab="Subset size",xlab="Iteration",pch=16,cex=0.5)
      text(x=1,y=1,i)
      plot(1:i,temp[1:i],ylim=c(0,1),xlim=c(0,animPlot),
           ylab="Temperature",xlab="Iteration",pch=16,cex=0.5)
      text(x=0,y=0,i)
      plot(1:i,score,ylim=c(-5,5),xlim=c(0,animPlot),
           ylab="Subset score",xlab="Iteration",pch=16,cex=0.5)
      abline(h=0,lty=2)
      text(x=0,y=-5,i)
      ani.pause()
    }
  }

  ### Return the results (subnetwork, size, score and p-value)
  Stat<-mean(workingTable[workingTable$state,]$score)
  subnetSize<-length(workingTable[workingTable$state,]$gene)
  pval<-1-pnorm(Stat,
                mean=nullDist[nullDist$k==length(activeNet),]$mu,
                sd=nullDist[nullDist$k==length(activeNet),]$sigma)
  cat(paste("\n\n  Subnetwork size:",subnetSize,
            "genes\n  Subnetwork score:",format(s,digits=4),
            "\n  p-value:",pval,"\n\n"))

  ret<-list(table=workingTable,score=s,size=subnetSize,pvalue=pval)
  invisible((ret))
}
