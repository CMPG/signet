#' Search for a high scoring subnetwork
#'
#' A simulated annealing algorithm to find the highest scoring
#' subnetwork in a graph.
#'
#' @param pathway A gene network, or a list of networks,
#'  in the \verb{graphNEL} format.
#' @param scores A data frame with two columns: gene identifiers list
#'  and associated scores.
#' @param nullDist A data frame with three columns : k, mu, sigma.
#' Can be obtained thanks to the \verb{nulldistribution()} function.
#' @param replicates Number of replicates per network.
#' @param iterations Number of iterations.
#' @param kmin The minimal size of a subnetwork.
#' @param verbose If \verb{TRUE}, displays text in the R console.
#'
#' @keywords subnetwork, simulated annealing, heuristics, search algorithm
#'
#' @return A list containing a table with genes, their state, their score;
#' the subnetwortk score and size and the p-value
#'
#' @export
#' @examples
#' require(signet)
#' data(keggPathways)
#' data(zScores)
#'
#' example1 <- searchSubnet(keggPathways[[1]], zScores)

searchSubnet<-function(pathway,
                       scores,
                       nullDist,
                       replicates = 1,
                       iterations = 5000,
                       subnetScore="ideker",
                       kmin = 2,
                       verbose = TRUE) {

  # check for packages =========================================================
  if (sum(installed.packages()[, 1]=="graph")==0) {
    stop("Package graph is required.")
  } else {
    requireNamespace("graph", quietly=TRUE)
  }

  #check for arguments =========================================================
  if (missing(pathway)) stop("Pathway data are missing.")
  if (missing(scores)) stop("Gene scores are missing.")
  if (missing(nullDist)) stop("Background distribution is missing.")
  colnames(scores) <- c("gene", "score")

  #check if list or unique pathway =============================================
  if (class(pathway) == "list") {
    all<-list()
    for (i in 1:length(pathway)) {
      print(names(pathway[i]))
      res<-try(
        searchSubnet(pathway[[i]],
                     scores=scores,
                     nullDist = nullDist,
                     iterations = iterations,
                     replicates = replicates,
                     verbose=verbose,
                     subnetScore=subnetScore,
                     kmin = kmin)
      )
      if (class(res)=="try-error") {
        res<-NA
      }
      all[[i]]<-res
      cat("\r  ", i, "/", length(pathway), " pathways analyzed.", sep="")
    }
    names(all)<-names(pathway)
    class(all)<-"signetList"
    invisible(all)

  } else if (class(pathway)!="graphNEL") {

    stop("Pathway is not a graphNEL object")

  } else {

    if (replicates>1) {

      if (verbose) {
        cat("Replicating: ")
      }
      allReturn<-replicate(replicates, {out<-searchSubnet(pathway=pathway,
                                                          scores = scores,
                                                          nullDist = nullDist,
                                                          replicates = -1,
                                                          iterations = iterations,
                                                          kmin = kmin,
                                                          verbose = FALSE,
                                                          subnetScore = subnetScore);if (verbose) cat("+");return(out)},simplify=FALSE)

      sz <- unlist(lapply(allReturn, function(x) return(x$subnet_size)))

      condS <- which(sz > 1) #plusieurs maxima possibles,
      allReturn <- allReturn[condS]

      vec <- unlist(lapply(allReturn, function(x) return(x$subnet_score)))

      if (length(vec)==0) {
        if (verbose) {
          cat("\nNo high-scoring subnetwork found\n\n")
        }
      } else {
        condV<-which(vec==max(vec))
        signetObject<-allReturn[condV][[1]]#on prend le premier si plusieurs
      }
    } else {

      # INITIALIZATION ==================================================
      signetObject<-createSignetObject(pathway,scores,iterations,5)

      if(length(nodes(signetObject$connected_comp))<5) {
        return(signetObject)
      }

      adjMatrix <- getAdjacencyMatrix(pathway=signetObject$connected_comp)
      boundaries <- NULL
      geneSampled <- sample(colnames(adjMatrix), 1)
      geneSampled <- c(geneSampled, sample(names(which(adjMatrix[, geneSampled]>0)))[1])

      # toggle state in final list
      signetObject$network[which(signetObject$network$gene%in%geneSampled), ]$active <- TRUE
      sumStat <- computeScore(signetObject, score = subnetScore)

      # SCORE COMPUTATION ===============================================
      s <- (sumStat-nullDist[nullDist$k == kmin, ]$mu)/nullDist[nullDist$k==kmin,]$sigma

      # OPTIMISATION ====================================================
      for (i in 1:iterations) {

        activeNet <- as.character(signetObject$network[signetObject$network$active,]$gene)
        adjSubgraph <- adjMatrix[activeNet,]

        if (length(dim(adjSubgraph))>0) {
          boundaries <- adjSubgraph[apply(adjSubgraph,1,sum)>0,]
          boundaries <- colnames(boundaries[,apply(boundaries,2,sum)>0])
          #remove boundaries already active
          boundaries <- boundaries[!boundaries %in% signetObject$network[signetObject$network$active,]$gene]
        } else {
          boundaries <- NULL
        }

        bla <- 0
        while (bla != 1) {
          ge <- sample(activeNet,1)
          test <- graph::subGraph(as.character(activeNet[activeNet != ge]),signetObject$connected_comp)

          bla <- length(graph::connComp(test))
          if (bla == 1) neighbours<-c(ge)
        }

        probneighbour <- 0
        sampNeighbour <- FALSE
        if (length(activeNet) > kmin) {
          probneighbour <- length(activeNet)/(length(boundaries)+length(activeNet))
          if (probneighbour > runif(1)) {
            sampNeighbour <- TRUE
          }
        }

        if (verbose) {
          if (length(activeNet) >= max(nullDist$k)) {
            stop(paste("No null distribution for subnetwork of size > ",max(nullDist$k)))
          }
        }

        if (length(activeNet) > 1 & length(unique(c(neighbours,boundaries))) > 0) {

          # Sample a new gene to toggle its state
          if (!sampNeighbour) newG <- as.character(sample(unique(c(boundaries)),1))
          else newG <- neighbours
          signetObject$network[which(signetObject$network$gene==newG),]$active <-
            !signetObject$network[which(signetObject$network$gene==newG),]$active

          # Compute the subnet score
          sumStat <- computeScore(signetObject, score = subnetScore)

          # Scale the subnet score
          s2 <- (sumStat-nullDist[nullDist$k == length(activeNet),]$mu)/nullDist[nullDist$k==length(activeNet),]$sigma

          # Keep or not the toggled gene, acceptance probability
          if (s2 < s) {
            prob <- exp((s2-s)/signetObject$simulated_annealing$temperature[i])
            if (prob > runif(1)) {
              s <- s2
              signetObject$simulated_annealing$score_evolution[i] <- s
              signetObject$simulated_annealing$size_evolution[i] <- sum(
                signetObject$network$active
              )
            } else {
              signetObject$network[which(signetObject$network$gene==newG),]$active <-
                !signetObject$network[which(signetObject$network$gene==newG),]$active
              signetObject$simulated_annealing$score_evolution[i] <- s
              signetObject$simulated_annealing$size_evolution[i] <- sum(
                signetObject$network$active
              )
            }
          } else if (s2 > s) {
            s <- s2
            signetObject$simulated_annealing$score_evolution[i] <- s
            signetObject$simulated_annealing$size_evolution[i] <- sum(
              signetObject$network$active
            )
          }
        }
      }

      ### Return the results (subnetwork, size, score)
      signetObject$subnet_size<-signetObject$simulated_annealing$size_evolution[iterations]
      signetObject$subnet_score<-signetObject$simulated_annealing$score_evolution[iterations]
      signetObject$aggregate_score<-sum(signetObject$network$score[signetObject$network$active])/sqrt(signetObject$subnet_size)
      signetObject$mean_score<-sum(signetObject$network$score[signetObject$network$active])/(signetObject$subnet_size)
      signetObject$subnet_genes<-signetObject$network$gene[signetObject$network$active]
    }
    invisible(signetObject)
  }
}

