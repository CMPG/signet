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
                       subnetScore="sum",
                       kmin = 2,
                       directed = FALSE,
                       verbose = TRUE,
                       burnin = 100,
                       animPlot = 0,
                       diagnostic = FALSE)
{

  # check for packages =========================================================
  if (sum(installed.packages()[,1]=="graph")==0) {
    stop("Package graph is not installed !")
  } else {
    requireNamespace("graph",quietly=TRUE)
  }

  #check for arguments =========================================================
  if (missing(pathway)) stop("Pathway is missing !")
  if (missing(scores)) stop("Gene scores are missing !")
  if (verbose) cat("  Running simulated annealing...\n")
  if (missing(nullDist)) {
    maximean <- TRUE
  }
  else maximean <- FALSE
  colnames(scores) <- c("gene","score")

  #check if list or unique pathway =============================================
  if (class(pathway) == "list") {
    all<-list()
    for (i in 1:length(pathway)) {
      res<-try(
        searchSubnet(pathway[[i]],
                     scores=scores,
                     nullDist = nullDist,
                     iterations = iterations,
                     replicates = replicates,
                     temperature = temperature,
                     diagnostic=diagnostic,
                     verbose=verbose)
      )
      if (class(res)=="try-error") {
        res<-NULL
      }
      all[[i]]<-res
      cat("\r  ",i,"/",length(pathway)," pathways analyzed.",sep="")
    }
    invisible(all)
  } else if (class(pathway)!="graphNEL") {
    stop("Pathway is not a graphNEL object")
  } else {
    if (replicates>1) {
      if (verbose) cat("  Replicating: ")
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
                                                         diagnostic=diagnostic);if (verbose) cat("+");return(out)},

                           simplify=FALSE)
      if (verbose) cat("\n  ... Done !")
      sz <- unlist(lapply(allReturn, function(x) return(x$size)))

      condS <- which(sz > 1) #plusieurs maxima possibles,
      allReturn <- allReturn[condS]

      vec <- unlist(lapply(allReturn, function(x) return(x$score)))

      if (length(vec)==0) {
        if (verbose) {
          cat("\n  No high-scoring subnetwork found\n\n")
        }
        ret <- NULL
      } else {
        condV<-which(vec==max(vec))
        ret<-allReturn[condV][[1]]#on prend le premier si plusieurs
      }
    } else {
      # remove connected components of size < threshold
      threshold<-10

      if (length(graph::nodes(pathway))==0) X<-NULL
      else  X<-unlist(graph::connComp(pathway)[!lapply(graph::connComp(pathway),
                                                       length)<threshold])
      if (is.null(X)) {
        if (verbose) {
          cat("\r  No connected component of size greater than 10 ")
        }
        ret<-NULL
      } else {
        pathway<-graph::subGraph(X,pathway)

        if (directed==FALSE) pathway<-graph::ugraph(pathway)

        names(scores)[[1]]<-"gene";names(scores)[[2]]<-"score"
        scores<-scores[scores$gene %in% graph::nodes(pathway),]

        if (maximean & verbose) {
          cat("  No background distribution provided. The mean will be maximized.\n")
        }

        # Initialization of the graph
        newpath<-graph::subGraph(as.character(scores[!is.na(scores$score),]$gene),pathway) #remove missing values

        X2<-unlist(graph::connComp(newpath)[!lapply(graph::connComp(newpath),length)<threshold])
        if (is.null(X2)) {
          if (verbose) cat("\r  No connected component of size greater than 10 ")
          ret <- NULL
        } else {
          adjMatrix <- getAdjacencyMatrix(pathway=newpath) #get the adjacency matrix of the graph

          #initialize the gene scores list (gene names, scores, and state)
          signetObject <- scores[!is.na(scores$score),]
          signetObject$state <- rep(FALSE,length(signetObject$score))

          temp <- temperatureFunction(iterations = iterations,
                                      param = temperature,
                                      burnin = burnin)

          if (animPlot > 0) {
            sizeEvolution <- array(NA,animPlot)
            scoreEvolution <- array(NA,animPlot)
          }
          boundaries<-NULL

          geneSampled<-sample(colnames(adjMatrix),1)
          geneSampled<-c(geneSampled,sample(names(which(adjMatrix[,geneSampled]>0)))[1])

          #toggle state in final list
          signetObject[which(signetObject$gene%in%geneSampled),]$state <- TRUE
          #     signetObject[which(!signetObject$gene%in%geneSampled),]$state <- FALSE
          sumStat <- computeScore(signetObject,score = subnetScore)
          if (subnetScore=="mean") sumStat<-mean(signetObject[signetObject$state,]$score)
          if (subnetScore=="sum") sumStat<-sum(signetObject[signetObject$state,]$score)
          if (subnetScore=="delta") sumStat<-mean(signetObject[signetObject$state,]$score)-mean(tail(sort(signetObject[!signetObject$state,]$score),5))

          # SCORE COMPUTATION ==================================================
          if (maximean) {
            s <- sumStat
          } else {
            s <- (sumStat-nullDist[nullDist$k==kmin,]$mu)/nullDist[nullDist$k==kmin,]$sigma
          }


          # OPTIMISATION ==================================================
          for (i in 1:iterations) {
            if (verbose & (100*i/iterations)%%10 ==0) {
              cat('\r  Searching for a high-scoring subnetwork...',paste(10*(100*i/iterations)%/%10,"%",sep =""))
              flush.console()
            }

            activeNet <- as.character(signetObject[signetObject$state,]$gene)
            adjSubgraph <- adjMatrix[activeNet,]

            if (length(dim(adjSubgraph))>0) {
              boundaries <- adjSubgraph[apply(adjSubgraph,1,sum)>0,]
              boundaries <- colnames(boundaries[,apply(boundaries,2,sum)>0])
              #remove boundaries already active
              boundaries <- boundaries[!boundaries %in% signetObject[signetObject$state,]$gene]
            } else boundaries <- NULL

            bla <- 0
            while (bla != 1) {
              ge <- sample(activeNet,1)
              test <- graph::subGraph(as.character(activeNet[activeNet != ge]),newpath)

              bla <- length(graph::connComp(test))
              if (bla == 1) neighbours<-c(ge)
            }

            probneighbour <- 0
            sampNeighbour <- FALSE
            if (length(activeNet) > 2) {
              probneighbour <- length(activeNet)/(length(boundaries)+length(activeNet))
              if (probneighbour > runif(1)) {
                sampNeighbour <- TRUE
              }
            }

            if (!maximean & verbose) {
              if (length(activeNet) >= max(nullDist$k)) stop(paste("No null distribution
                                   for subnetwork
                                   of size > ",max(nullDist$k)))
            }

            if (length(activeNet) > 1 & length(unique(c(neighbours,boundaries))) > 0) {
              #Sample a new gene to toggle its state
              if (!sampNeighbour) newG <- as.character(sample(unique(c(boundaries)),1))
              else newG <- neighbours
              signetObject[which(signetObject$gene==newG),]$state <-
                !signetObject[which(signetObject$gene==newG),]$state

              # Compute the subnet score
              if (subnetScore=="mean") {
                sumStat <- mean(signetObject[signetObject$state,]$score)
              } else if (subnetScore=="sum") {
                sumStat <- sum(signetObject[signetObject$state,]$score)
              } else if (subnetScore=="delta") {
                sumStat <- mean(signetObject[signetObject$state,]$score)-mean(tail(sort(signetObject[!signetObject$state,]$score),5))
              }

              # Scale the subnet score
              if (maximean) {
                s2 <- sumStat
              } else {
                s2 <- (sumStat-nullDist[nullDist$k == length(activeNet),]$mu)/nullDist[nullDist$k==length(activeNet),]$sigma
              }

              # Keep or not the toggled gene, acceptance probability
              if (s2 < s) {
                prob <- exp((s2-s)/temp[i])
                if (prob > runif(1)) {
                  s <- s2
                } else {
                  signetObject[which(signetObject$gene==newG),]$state <- !signetObject[which(signetObject$gene==newG),]$state
                }
              } else if (s2 > s) {
                s <- s2
              }
            }

            if (diagnostic | (animPlot > 0 & i <= animPlot)) {
              sizeEvolution[i] <- length(signetObject[signetObject$state,]$gene)
              scoreEvolution[i] <- s

              if (animPlot > 0 & i <= animPlot & i %% 100 == 0) {
                m <- rbind(c(0,1,1,0),c(2,2,3,3))
                layout(m)

                plotSubnet(newpath,signetObject[signetObject$state,]$gene)

                par(mar=c(5,5,5,5))
                plot(1:i,
                     sizeEvolution[!is.na(sizeEvolution)],
                     ylim=c(1,length(graph::nodes(newpath))),
                     xlim=c(1,animPlot),
                     ylab="Subnetwork size",
                     xlab="Iteration",
                     pch=16,
                     cex=0.5,
                     col="firebrick",
                     type="l")

                plot(1:i,scoreEvolution[!is.na(scoreEvolution)],
                     ylim=c(-5,15),
                     xlim=c(0,animPlot),
                     ylab="Subnetwork score",
                     xlab="Iteration",
                     pch=16,cex=0.5,
                     col="dodgerblue",
                     type="l")

                requireNamespace("animation",quietly=TRUE)
                animation::ani.pause()
              }
            }
            if (diagnostic) {
              if (i==1) size <- NULL
              size<-c(size,length(signetObject[signetObject$state,]$gene))
              if (i==1) score <- NULL
              score <- c(score,s)
              if (i==iterations) {
                par(mfrow = c(1,2))

                x <- 1:iterations
                y <- sizeEvolution
                z <- scoreEvolution
                par(mar = c(5, 5, 5, 5))  # Leave space for z axis

                plot(x, z, type = "l",lwd = 1,cex = 0.2,pch = 16,
                     col = "firebrick",ylab = "Subnetwork score",
                     xlab = "Iterations") # first plot

                par(new = TRUE)
                plot(1:iterations, temp, type = "l", lwd=1,
                     xlab = "", ylab = "", axes=F, lty=2,
                     col = "grey")
                axis(side = 4, at = pretty(range(temp)))
                mtext("Temperature", side = 4, line = 3)

                plot(x, y, type = "l", cex = 0.5,ylim=c(0,max(y)),
                     xlab = "Iterations", ylab = "Subnetwork size",
                     pch=16,lwd=1,col="dodgerblue")

                par(new=TRUE)
                plot(1:iterations,temp, type="l",lwd=1,
                     xlab="", ylab="",axes=F,lty=2,
                     col="grey")
                axis(side=4, at = pretty(range(temp)))
                mtext("Temperature", side=4, line=3)

                par(mfrow = c(1,1))
              }
            }
          }

          ### Return the results (subnetwork, size, score and p-value)
          Stat<-mean(signetObject[signetObject$state,]$score)
          subnetSize<-length(signetObject[signetObject$state,]$gene)
          ret<-list(table=signetObject,score=s,size=subnetSize)
        }
      }
    }
    if (verbose & !is.null(ret)) {
      cat(paste("\n\n  Subnetwork size:",ret$size,
                "genes\n  Subnetwork score:",format(ret$score,digits=4),
                "\n"))
    }
    invisible(ret)
  }
}
