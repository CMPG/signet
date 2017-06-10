#' Search for a high scoring subnetwork
#'
#' A simulated annealing algorithm to find the highest scoring
#' subnetwork in a graph.
#'
#' @param pathway A gene network, or a list of networks,
#'  in the \verb{graphNEL} format.
#' @param scores A data frame with two columns: gene identifiers list
#'  and associated scores.
#' @param background A data frame with three columns : k, mu, sigma.
#' @param iterations Number of iterations.
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

searchSubnet <- function(pathway,
                         scores,
                         background,
                         iterations = 5000,
                         subnetScore="ideker",
                         verbose = TRUE) {

  # check for packages =========================================================
  if (sum(installed.packages()[, 1]=="graph")==0) {
    stop("Package graph is required.")
  } else {
    requireNamespace("graph", quietly=TRUE)
    requireNamespace("RBGL", quietly=TRUE)
  }

  #check for arguments =========================================================
  if (missing(pathway)) stop("Pathway data are missing.")
  if (missing(scores)) stop("Gene scores are missing.")
  colnames(scores) <- c("gene", "score")

  if (class(pathway) == "list") {

    uniqGenes<-unique(unlist(sapply(pathway,function(x) unique(nodes(x)))))
    scores<-scores[scores$gene %in% uniqGenes,]

  } else if (class(pathway)!="graphNEL") {

    uniqGenes<-unique(unlist(nodes(pathway)))
    scores<-scores[scores$gene %in% uniqGenes,]

  }

  if (missing(background)) {
    background <- data.frame(k=1:100,
                             mu=sqrt(1:100)*mean(scores$score,na.rm=TRUE),
                             sigma=rep(sd(scores$score,na.rm=TRUE),100))
  }

  #check if list or unique pathway =============================================
  if (class(pathway) == "list") {
    all<-list()
    for (i in 1:length(pathway)) {
      cat(paste("\n  Analyzing pathway: ",names(pathway[i])))
      res<-try(
        searchSubnet(pathway[[i]],
                     scores=scores,
                     background = background,
                     iterations = iterations,
                     verbose=FALSE,
                     subnetScore=subnetScore)
      )
      if (class(res)=="try-error") {
        res<-NA
      }
      all[[i]]<-res
      cat("\n  ", i, "/", length(pathway), " pathways analyzed.", sep="")
    }
    names(all)<-names(pathway)
    class(all)<-"signetList"
    invisible(all)

  } else if (class(pathway)!="graphNEL") {

    stop("Pathway is not a graphNEL object")

  } else {

    # INITIALIZATION ==================================================
    sigObj<-createSignetObject(pathway,scores,iterations,5)

    if(length(graph::nodes(sigObj$connected_comp))<5) {
      return(sigObj)
    }

    adjMatrix <- getAdjacencyMatrix(pathway=sigObj$connected_comp)

    # sample two connected genes
    geneSampled <- sample(colnames(adjMatrix), 1)
    geneSampled <- c(geneSampled, sample(names(which(adjMatrix[, geneSampled]>0)))[1])

    # toggle state of the two starting genes
    sigObj$network[which(sigObj$network$gene%in%geneSampled), ]$active <- TRUE

    # SCORE COMPUTATION ===============================================
    sumStat <- computeScore(sigObj, score = subnetScore)
    s <- (sumStat-background[background$k == 2, ]$mu)/background[background$k==2,]$sigma

    if(verbose) cat("  Running simulated annealing...\n")
    rn <- runif(iterations)
    rn2 <- runif(iterations)

    boundaries <- NULL

    # OPTIMISATION ====================================================
    for (i in 1:iterations) {

      if(verbose & i %% 100 == 0) cat(paste("\r  Iteration ",i))

      activeNet <- sigObj$network[sigObj$network$active,]$gene
      activeNet <- as.character(activeNet)
      adjSubgraph <- adjMatrix[activeNet,]

      if (length(dim(adjSubgraph))>0) {

        boundaries <- adjSubgraph[apply(adjSubgraph,1,sum)>0,]
        boundaries <- colnames(boundaries[,apply(boundaries,2,sum)>0])
        # remove boundaries already active
        boundaries <- boundaries[!boundaries %in% activeNet]

      } else {

        boundaries <- NULL

      }

      neighbours <- NULL
      probneighbour <- 0
      sampNeighbour <- FALSE

      if (length(activeNet) > 2) {

        probneighbour <- 0.5

        if (probneighbour > rn2[i]) {

          sampNeighbour <- TRUE

          if (length(activeNet) <= 10) {

            CComp <- list(1,2)

            while (length(CComp) != 1) { #bottleneck

              ge <- sample(activeNet,1)
              test <- graph::subGraph(as.character(activeNet[activeNet != ge]),
                                      sigObj$connected_comp)
              CComp <- RBGL::connectedComp(test)

              if (length(CComp) == 1) {
                neighbours<-c(ge)
              }

            }

          } else {

            ge <- sample(activeNet,1)
            test <- graph::subGraph(as.character(activeNet[activeNet != ge]),
                                    sigObj$connected_comp)
            CComp <- RBGL::connectedComp(test)

            if(max(sapply(CComp,length)) <= 10) {
              sampNeighbour <- FALSE
            }

            neighbours<-c(ge)#added

          }
        }
      }

      if (verbose & length(activeNet) >= max(background$k)) {
        stop(paste("No background distribution for size > ",max(background$k)))
      }

      if (length(activeNet) > 1 & length(unique(c(neighbours,boundaries))) > 0) {

        # Sample a new gene to toggle its state
        if (!sampNeighbour) newG <- as.character(sample(unique(c(boundaries)),1))
        else newG <- neighbours

        sigObj$network[which(sigObj$network$gene%in%newG),]$active <-
          !sigObj$network[which(sigObj$network$gene%in%newG),]$active


        # Compute the subnet score
        if (!sampNeighbour){

          sumStat <- computeScore(sigObj, score = subnetScore)

          # Scale the subnet score
          s2 <- (sumStat-background[background$k == length(activeNet),]$mu)/
            background[background$k==length(activeNet),]$sigma

        } else {

          if(length(CComp)==1){

            sc2<-computeScore(sigObj, score = subnetScore)
            cK<-sum(sigObj$network$active)
            s2<-(sc2-background[background$k==cK,]$mu)/
              (background[background$k==cK,]$sigma)

          } else {

            scs<-lapply(CComp,function(x){
              cK<-length(x)

              if(cK >= 5) {
                scn<-(1/(sqrt(cK)))*sum(sigObj$network$score[
                  as.character(sigObj$network$gene)%in%x]
                  )
                scn<-(scn-background[background$k==cK,]$mu)/
                  (background[background$k==cK,]$sigma)
              } else {
                scn <- NA
              }

            }) #compute score for each CC

            s2<-max(unlist(scs),na.rm=TRUE)
            GL<-CComp[[which(unlist(scs)==s2)[1]]]
            sigObj$network[!(as.character(sigObj$network$gene) %in% GL),]$active <- FALSE

          }
        }

        # Keep or not the toggled gene, acceptance probability
        if (s2 < s) {

          #acceptance probability
          prob <- exp((s2-s)/sigObj$simulated_annealing$temperature[i])

          if (prob > rn[i]) {
            s <- s2
            sigObj$simulated_annealing$score_evolution[i] <- s
            sigObj$simulated_annealing$size_evolution[i] <- sum(
              sigObj$network$active
            )
          } else {
            sigObj$network[which(sigObj$network$gene%in%newG),]$active <-
              !sigObj$network[which(sigObj$network$gene%in%newG),]$active
            sigObj$simulated_annealing$score_evolution[i] <- s
            sigObj$simulated_annealing$size_evolution[i] <- sum(
              sigObj$network$active
            )
          }
        } else if (s2 >= s) {
          s <- s2
          sigObj$simulated_annealing$score_evolution[i] <- s
          sigObj$simulated_annealing$size_evolution[i] <- sum(
            sigObj$network$active
          )
        }
      }
    }

    ### Return the results (subnetwork, size, score)
    sigObj$subnet_size<-sigObj$simulated_annealing$size_evolution[iterations]
    sigObj$subnet_score<-sigObj$simulated_annealing$score_evolution[iterations]
    sigObj$aggregate_score<-sum(sigObj$network$score[sigObj$network$active])/sqrt(sigObj$subnet_size)
    sigObj$mean_score<-sum(sigObj$network$score[sigObj$network$active])/(sigObj$subnet_size)
    sigObj$subnet_genes<-sigObj$network$gene[sigObj$network$active]
    invisible(sigObj)
  }
}
