#' Search for a high scoring subnetwork
#'
#' A simulated annealing algorithm is used to research a high scoring
#' subnetwork.Implements also a test id a null distribution is provided.
#'
#' @param pathway A gene network, or a list of networks,
#'  in the \verb{graphNEL} format.
#' @param scores A data frame with two columns: gene identifiers list
#'  and associated scores.
#' @param nullDist A data frame with three columns : k, mu, sigma.
#' Can be obtained thanks to the \verb{nulldistribution()} function.
#' @param replicates Number of replicates per network.
#' @param iterations Number of iterations.
#' @param temperature The temperature function parameter.
#' @param kmin The minimal size of a subnetwork.
#' @param directed If \verb{TRUE}, considers the edges direction, i.e. cannot go back.
#' @param verbose If \verb{TRUE}, displays text in the R console.
#' @param burnin Number of iterations before the temperature begins to decrease.
#' @param animPlot For development purposes.
#' @param diagnostic If TRUE, plots the evolution of different stats along the run.
#'
#' @keywords subnetwork, simulated annealing, heuristics, search algorithm
#'
#' @return A list containing a table with genes, their state, their score;
#' the subnetwortk score and size and the p-value
#'
#' @export
#'
#' @examples
#' require(signet)
#' data(keggPathways)
#' data(zScores)
#'
#' example1 <- searchSubnet(keggPathways[[1]],zScores)
#'
#' #data() #a previously generated null distribution
#' #example2 <- searchSubnet(keggPathways[[1]],zScores,null)

searchSubnet<-function(pathway,
                       scores,
                       nullDist,
                       replicates = 1,
                       iterations = 2000,
                       temperature = 0.995,
                       kmin = 3,
                       directed = FALSE,
                       verbose = TRUE,
                       burnin = 100,
                       animPlot = 0,
                       diagnostic=FALSE)
{
  #check for the package graph and load it
  if(sum(installed.packages()[,1]=="graph")==0)
  {
    stop("Package graph is not installed !")
  }
  else
  {
    requireNamespace("graph",quietly=TRUE)
  }
  #check for arguments
  if(missing(pathway)) stop("Pathway is missing !")
  if(missing(scores)) stop("Gene scores are missing !")
  if(verbose) cat("  Running simulated annealing...\n")
  if(missing(nullDist))
  {
    maximean<-TRUE
  }
  else maximean<-FALSE
  colnames(scores)<-c("gene","score")

  # if(diagnostic & replicates>0) par(mfrow=c(replicates,2))

  if(replicates>1)
  {
    if(verbose) cat("  Replicating: ")
    allReturn<-replicate(replicates,{out<-searchSubnet(pathway=pathway,
                      scores = scores,
                      nullDist = nullDist,
                      replicates = -1,
                      iterations = iterations,
                      temperature = temperature,
                      kmin = kmin,
                      directed = directed,
                      verbose = FALSE,
                      burnin = burnin,
                      animPlot = 0,
                      diagnostic=diagnostic);cat("+");return(out)},

                      simplify=FALSE)
    cat("\n  ... Done !")
    sz<-unlist(lapply(allReturn,function(x) return(x$size)))


    condS<-which(sz > kmin) #plusieurs maxima possibles,
    allReturn<-allReturn[condS]

    vec<-unlist(lapply(allReturn,function(x) return(x$score)))

    if(length(vec)==0){
      cat("\n  No high-scoring subnetwork found\n\n")
      ret <- NULL
    }
    else
    {
      condV<-which(vec==max(vec))
      ret<-allReturn[condV][[1]]#on prend le premier si plusieurs
    }
  }
  else
  {
  #remove conn components of size < threshold
  threshold<-10
  X<-unlist(graph::connComp(pathway)[!lapply(graph::connComp(pathway),length)<threshold])
  if(is.null(X))
  {
    cat("\r  No connected component of size greater than 10 ")
    ret<-NULL
  }
  else
  {
    pathway<-graph::subGraph(X,pathway)

    if(directed==FALSE) pathway<-graph::ugraph(pathway)

    names(scores)[[1]]<-"gene";names(scores)[[2]]<-"score"
    scores<-scores[scores$gene %in% graph::nodes(pathway),]

    if(maximean & verbose)
    {
      cat("  No null distribution provided. The mean will be maximized.\n")
    }

    ### Initialization of the graph

    newpath<-graph::subGraph(as.character(scores[!is.na(scores$score),]$gene),pathway) #remove missing values
    adjMatrix<-getAdjacencyMatrix(pathway=newpath) #get the adjacency matrix of the graph

    #initialize the gene scores list (gene names, scores, and state)
    workingTable <- scores[!is.na(scores$score),]
    workingTable$state <- rep(FALSE,length(workingTable$score))

    #initialize the temperature function
    temp<-temperatureFunction(iterations = iterations, param = temperature, burnin = burnin)

    #c. select the initial subset of size kmin
    boundaries<-NULL
    geneSampled<-NULL
    for(ii in 1:kmin)
    {
      bu<-0
      countW<-0
      gene<-NULL
      while(sum(bu)==0)
      {
        if(ii == 1){ gene<-sample(colnames(adjMatrix),1)}## 1. 1 gene au hasard
        else{ gene<-sample(boundaries,1)}
        bu<-adjMatrix[gene,]
        countW<-countW+1
        if(countW==100){ret<-NULL;break()}
      }
      if(countW==100){break()}
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

    if(countW<100){

    #toggle state in final list
    workingTable[which(workingTable$gene%in%geneSampled),]$state <- TRUE
    workingTable[which(!workingTable$gene%in%geneSampled),]$state <- FALSE
    sumStat<-mean(workingTable[workingTable$state,]$score)

    if(maximean)
    {
      s <- sumStat
    }
    else
    {
      s<-(sumStat-nullDist[nullDist$k==kmin,]$mu)/nullDist[nullDist$k==kmin,]$sigma
    }


    ###2. OPTIMISATION
    ######

    for(i in 1:iterations)
    {

      if(verbose & (100*i/iterations)%%10==0)#displays the progession every 10%
      {
        cat('\r  Searching for a high-scoring subnetwork...',paste(10*(100*i/iterations)%/%10,"%",sep=""))
        flush.console()
      }

      ###
      activeNet<-as.character(workingTable[workingTable$state,]$gene)
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
        test<-graph::subGraph(as.character(activeNet[activeNet!=ge]),newpath)

        bla<-length(graph::connComp(test))
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

      if(!maximean & verbose)
      {
        if(length(activeNet)>=max(nullDist$k)) stop(paste("No null distribution for subnetwork
                                   of size > ",max(nullDist$k)))
      }

      if(length(activeNet)>=kmin & length(unique(c(neighbours,boundaries)))>0)
      {
        #Sample a new gene to toggle its state
        if(sampBound) newG<-as.character(sample(unique(c(neighbours,boundaries)),1))
        else newG<-neighbours
        workingTable[which(workingTable$gene==newG),]$state <- !workingTable[which(workingTable$gene==newG),]$state

        #Compute the subnet score
        sumStat<-mean(workingTable[workingTable$state,]$score)

        if(maximean)
        {
          s2 <- sumStat
        }
        else  s2<-(sumStat-nullDist[nullDist$k==length(activeNet),]$mu)/nullDist[nullDist$k==length(activeNet),]$sigma

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
             ylab="Subnetwork size",xlab="Iteration",pch=16,cex=0.5)
        text(x=1,y=1,i)
        plot(1:i,temp[1:i],ylim=c(0,1),xlim=c(0,animPlot),
             ylab="Temperature",xlab="Iteration",pch=16,cex=0.5)
        text(x=0,y=0,i)
        plot(1:i,score,ylim=c(-5,5),xlim=c(0,animPlot),
             ylab="Subnetwork score",xlab="Iteration",pch=16,cex=0.5)
        abline(h=0,lty=2)
        text(x=0,y=-5,i)
        requireNamespace("animation",quietly=TRUE)
        animation::ani.pause()
      }
      if(diagnostic)
      {
        if(i==1)size<-NULL
        size<-c(size,length(workingTable[workingTable$state,]$gene))
        if(i==1)score<-NULL
        score<-c(score,s)
        if(i==iterations)
        {

        # par(mfrow=c(2,2))

        par(mfrow=c(1,1))

        x <- 1:iterations
        y <- size

        z <- score
        par(mar = c(4, 4, 4, 4))  # Leave space for z axis

        plot(x, z, type="p",pch=16,cex=0.6,
             col="red",ylab="Subnetwork score",
             xlab="Iterations") # first plot

        par(new = TRUE)
        plot(x, y, type = "l", axes = FALSE,cex=0.8,
             bty = "n", xlab = "", ylab = "",
             pch=16,lwd=1,col="blue")
        axis(side=4, at = pretty(range(y)))
        mtext("Subnetwork size", side=4, line=3)

        par(new=TRUE)

        plot(1:iterations,temp, type="l",lwd=1,
             xlab="", ylab="",axes=F,
             col="grey")
#         axis(1, labels=NA,at=c(0,5))
#         axis(2, labels=NA,at=c(0,150))


#         plot(1:iterations,size,ylim=c(0,max(size)),xlim=c(0,iterations),
#              ylab="Subnetwork size",xlab="Iteration",pch=16,cex=0.5)
#         plot(1:iterations,temp,ylim=c(0,1),xlim=c(0,iterations),
#              ylab="Temperature",xlab="Iteration",pch=16,cex=0.5)
#         plot(1:iterations,score,ylim=c(min(score),max(score)),xlim=c(0,iterations),
#              ylab="Subnetwork score",xlab="Iteration",pch=16,cex=0.5)
#         abline(h=0,lty=2)

#         plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n',
#              xaxt = 'n', yaxt = 'n')
#
#         text(x = 0.5, y = 0.5, paste("Subnet size : ",size[iterations],"\n",
#                              "Score :",format(score[iterations],digits=4)),
#         cex = 1.5, col = "darkred")


        }
      }
    }

    ### Return the results (subnetwork, size, score and p-value)
    Stat<-mean(workingTable[workingTable$state,]$score)
    subnetSize<-length(workingTable[workingTable$state,]$gene)
    if(!missing(nullDist))
    {
    pval<-1-pnorm(Stat,
                  mean=nullDist[nullDist$k==length(activeNet),]$mu,
                  sd=nullDist[nullDist$k==length(activeNet),]$sigma)
    }
    else pval <- NA
    ret<-list(table=workingTable,score=s,size=subnetSize,pvalue=pval)

    }
  }
  }
  if(verbose & !is.null(ret))
  {
    cat(paste("\n\n  Subnetwork size:",ret$size,
              "genes\n  Subnetwork score:",format(ret$score,digits=4),
              "\n  p-value:",ret$pvalue,"\n\n"))
  }

  invisible(ret)

}
