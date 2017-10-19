#' Search for a high scoring subnetwork
#'
#' A simulated annealing algorithm to find the highest scoring
#' subnetwork within a graph.
#'
#' @param pathway A gene network, or a list of gene networks,
#'  in the \verb{graphNEL} format.
#' @param scores A data frame with two columns: gene identifiers list (IDs have
#' to be the same as for the pathways, e.g. Entrez) and associated scores.
#' @param background Optional. A data frame generated using the
#' \verb{backgroundDist} function. Is missing, the
#' theoretical expectation is computed (see Gouy et al., 2017).
#' @param iterations Number of iterations.
#' @param verbose If \verb{TRUE}, displays text in the R console.
#'
#' @keywords subnetwork, simulated annealing, heuristics, search algorithm
#'
#' @return A \verb{signet} object or a list of \verb{signet} objects. Each
#' \verb{signet} object consists in a table with gene IDs, their state,
#' their score; the subnetwork score and size and the p-value.
#'
#' @export
#' @examples
#' data(daub2013) #pathways and gene scores from Daub et al. (2013).
#'
#' #run the search in the first two pathways with 1000 iterations
#' example <- searchSubnet(kegg_human[1:2], scores, iterations = 1000)
#' \dontrun{
#'
#' #run the search in all the pathways with 5000 iterations (default)
#' example <- searchSubnet(kegg_human, scores)
#' }

searchSubnet <- function(pathway,
                         scores,
                         background,
                         iterations = 5000,
                         verbose = TRUE) {

    # check for packages =======================================================

    if (sum(installed.packages()[, 1]=="graph")==0) {
        stop("Package graph is required.")
    } else {
        requireNamespace("graph", quietly=TRUE)
        requireNamespace("RBGL", quietly=TRUE)
        requireNamespace("stats", quietly=TRUE)
    }

    #check for arguments =======================================================

    if (missing(pathway)) stop("Pathway data are missing.")
    if (missing(scores)) stop("Gene scores are missing.")
    colnames(scores) <- c("gene", "score")

    if (missing(background)) {
        background <- data.frame(k = 1:200,
                                 mu = sqrt(1:200)*mean(scores$score,na.rm=TRUE),
                                 sigma = rep(sd(scores$score,na.rm=TRUE), 200))
    }

    if (class(pathway) == "list" & length(pathway) > 100) {

        uniqGenes <- unique(unlist(sapply(pathway,function(x) {
            unique(graph::nodes(x))
        })))
        scores <- scores[scores$gene %in% uniqGenes,]

    }

    #check if list or unique pathway ===========================================

    if (class(pathway) == "list") {
        all <- list()
        for (i in 1:length(pathway)) {
            message(paste("\r\r  Analyzing pathway ",
                          i,
                          "/",
                          length(pathway),
                          sep=""),
                    appendLF = TRUE)

            res <- try(
                searchSubnet(pathway[[i]],
                             scores=scores,
                             background = background,
                             iterations = iterations,
                             verbose=FALSE),
                silent = TRUE
            )

            if (class(res)=="try-error") {
                res <- NA
            }

            all[[i]] <- res

        }
        names(all) <- names(pathway)
        class(all) <- "signetList"

        if(verbose) message("  Done!")

        invisible(all)

    } else if (class(pathway)!="graphNEL") {

        stop("Pathway is not a graphNEL object")

    } else {

        # INITIALIZATION =======================================================

        sigObj <- createSignetObject(pathway,scores,iterations)

        uniqGenes <- unique(unlist(graph::nodes(pathway)))
        scores <- scores[scores$gene %in% uniqGenes,]

        if(length(graph::nodes(sigObj$connected_comp))<5) {
            return(sigObj)
        }

        adjMatrix <- getAdjacencyMatrix(pathway=sigObj$connected_comp)

        # sample two connected genes
        geneSampled <- sample(colnames(adjMatrix), 1)
        geneSampled <- c(geneSampled,
                         sample(names(which(adjMatrix[, geneSampled]>0)))[1])

        # toggle state of the two starting genes
        sigObj$network[which(sigObj$network$gene%in%geneSampled), ]$active <-
            TRUE

        # SCORE COMPUTATION ====================================================

        sumStat <- computeScore(sigObj, score = "ideker")
        s <- (sumStat-background[background$k == 2, ]$mu)/
            background[background$k==2,]$sigma

        if(verbose) message("  Running simulated annealing...\n",
                            appendLF = FALSE)
        rn <- runif(iterations)
        rn2 <- runif(iterations)

        boundaries <- NULL

        # OPTIMISATION =========================================================

        for (i in 1:iterations) {

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
            lenActive <- length(activeNet)

            if (lenActive > 2) {

                probneighbour <- 0.5

                if (probneighbour > rn2[i]) {

                    sampNeighbour <- TRUE

                    if (lenActive <= 10) {

                        CComp <- list(1,2)

                        while (length(CComp) != 1) { #bottleneck

                            ge <- sample(activeNet,1)
                            test <- graph::subGraph(as.character(
                                activeNet[activeNet != ge]),
                                                    sigObj$connected_comp)
                            CComp <- RBGL::connectedComp(test)

                            if (length(CComp) == 1) {
                                neighbours <- c(ge)
                            }

                        }

                    } else {

                        ge <- sample(activeNet,1)
                        test <- graph::subGraph(as.character(
                            activeNet[activeNet != ge]),
                                                sigObj$connected_comp)
                        CComp <- RBGL::connectedComp(test)

                        if(max(sapply(CComp,length)) < 5) {
                            sampNeighbour <- FALSE
                        }

                        neighbours <- c(ge)

                    }
                }
            }

            if (verbose & lenActive >= max(background$k)) {

                stop(paste("No background distribution for size > ",
                           max(background$k)))

            }

            if (lenActive > 1 & length(unique(c(neighbours,boundaries))) > 0) {

                # Sample a new gene to toggle its state
                if (!sampNeighbour) {
                    newG <- as.character(sample(unique(c(boundaries)),1))
                } else {
                    newG <- neighbours
                }

                sigObj$network[which(sigObj$network$gene%in%newG),]$active  <-
                    !sigObj$network[which(sigObj$network$gene%in%newG),]$active

                # Compute the subnet score
                if (!sampNeighbour) {

                    sumStat <- computeScore(sigObj, score = "ideker")

                    # Scale the subnet score
                    lenActive <- sum(sigObj$network$active)
                    s2 <- (sumStat-background[background$k == lenActive,]$mu)/
                        background[background$k==lenActive,]$sigma

                } else {

                    if(length(CComp)==1){

                        sc2 <- computeScore(sigObj, score = "ideker")
                        cK <- sum(sigObj$network$active)
                        s2 <- (sc2-background[background$k==cK,]$mu)/
                            (background[background$k==cK,]$sigma)

                    } else {

                        scs <- lapply(CComp,function(x){
                            cK <- length(x)

                            if(cK >= 5) {
                                scn <- (1/(sqrt(cK)))*sum(sigObj$network$score[
                                    as.character(sigObj$network$gene)%in%x]
                                )
                                scn <- (scn-background[background$k==cK,]$mu)/
                                    (background[background$k==cK,]$sigma)
                            } else {
                                scn <- NA
                            }

                        }) #compute score for each CC

                        s2 <- max(unlist(scs),na.rm=TRUE)
                        GL <- CComp[[which(unlist(scs)==s2)[1]]]
                        nG <- as.character(sigObj$network$gene)
                        sigObj$network[!(nG %in% GL), ]$active <- FALSE

                    }
                }

                # Keep or not the toggled gene, acceptance probability
                if (s2 < s) {

                    #acceptance probability
                    prob <- exp((s2-s)/
                                    sigObj$simulated_annealing$temperature[i])
                    if(i > iterations-iterations/10) prob <- 0

                    if (prob > rn[i]) {
                        s <- s2
                        sigObj$simulated_annealing$score_evolution[i] <- s
                        sigObj$simulated_annealing$size_evolution[i] <- sum(
                            sigObj$network$active
                        )
                    } else {
                        cond <- which(sigObj$network$gene %in% newG)

                        sigObj$network[cond, ]$active  <-
                            !sigObj$network[cond,]$active

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
        sigObj$subnet_size <-
            sigObj$simulated_annealing$size_evolution[iterations]

        sigObj$subnet_score <-
            sigObj$simulated_annealing$score_evolution[iterations]

        sigObj$aggregate_score <-
            sum(sigObj$network$score[sigObj$network$active])/
            sqrt(sigObj$subnet_size)

        sigObj$mean_score <-
            sum(sigObj$network$score[sigObj$network$active])/
            (sigObj$subnet_size)

        sigObj$subnet_genes <- sigObj$network$gene[sigObj$network$active]

        invisible(sigObj)
    }
}
