#' Search for a high scoring subnetwork
#'
#' A simulated annealing algorithm to find the highest scoring
#' subnetwork within a graph.
#'
#' @param pathway A gene network, or a list of gene networks,
#'  in the \verb{graphNEL} format.
#' @param scores A data frame with two columns: gene identifiers list (IDs have
#' to be the same as for the pathways, e.g. Entrez) and associated scores.
#' @param iterations Number of iterations.
#' @param background For development purposes.
#'
#' @keywords subnetwork, simulated annealing, heuristics, search algorithm
#'
#' @return A \verb{signet} object or a list of \verb{signet} objects. Each
#' \verb{signet} object consists in a table with gene IDs, their state,
#' their score; the subnetwork score and size and the p-value.
#'
#' @export
#' @examples
#'
#' # Get KEGG pathways from the package graphite:
#' library(graphite)
#' kegg <- pathways("hsapiens", "kegg")
#' kegg_human <- lapply(kegg[1:5], pathwayGraph)
#'
#' data(daub2013) # load the gene scores from Daub et al. (2013)
#'
#' #run the search in all the pathways with 2500 iterations (default)
#' example <- searchSubnet(kegg_human, scores)
#' summary(example)

searchSubnet <- function(pathway, scores, iterations = 1000, background) {

    # check for packages =======================================================

    if (sum(installed.packages()[, 1] == "graph") == 0) {
        stop("Package graph is required.")
    } else {
        requireNamespace("graph", quietly = TRUE)
        requireNamespace("RBGL", quietly = TRUE)
        requireNamespace("stats", quietly = TRUE)
    }

    #check for arguments =======================================================

    if (missing(pathway)) stop("Pathway data are missing.")
    if (missing(scores)) stop("Gene scores are missing.")
    colnames(scores) <- c("gene", "score")

    if (class(pathway) == "list" & length(pathway) > 100) {

        uniqGenes <- unique(unlist(sapply(pathway, function(x) {
            unique(graph::nodes(x))
        })))
        scores <- scores[scores$gene %in% uniqGenes, ]

    }

    if (missing(background)) {
        sizeP <- 200
        background <- data.frame(
            mu = sqrt(seq_len(sizeP)) * mean(scores$score, na.rm = TRUE),
            sigma = rep(sd(scores$score, na.rm = TRUE), sizeP)
        )
    }
    bk <- background

    #check if list or unique pathway ===========================================

    if (class(pathway) == "list") {

        message("  Running simulated annealing...")
        all <- lapply_pb(pathway, function(x) {

            res <- try(
                searchSubnet(
                    x,
                    scores = scores,
                    background = bk,
                    iterations = iterations
                ),
                silent = TRUE
            )
            if (class(res) == "try-error") res <- NA
            res

        })

        cond <- lapply(all, function(x) class(x)=="Signet")
        all <- .SignetList(all[unlist(cond)])

        invisible(all)

    } else if (class(pathway) != "graphNEL") {

        stop("Pathway is not a graphNEL object")

    } else {

        # INITIALIZATION =======================================================
        sigObj <- .Signet(pathway, scores, iterations)

        uniqGenes <- unique(unlist(graph::nodes(pathway)))
        scores <- scores[scores$gene %in% uniqGenes, ]

        if(length(graph::nodes(sigObj@connected_comp)) < 5) return(sigObj)

        adjMatrix <- getAdjacencyMatrix(pathway = sigObj@connected_comp)

        # sample two connected genes
        geneSampled <- sample(colnames(adjMatrix), 1)
        geneSampled <- c(
            geneSampled,
            sample(names(which(adjMatrix[, geneSampled] > 0)))[1]
        )

        # toggle state of the two starting genes
        sigObj@network[
            which(sigObj@network$gene%in%geneSampled),
            ]$active <- TRUE

        # SCORE COMPUTATION ====================================================

        sumStat <- computeScore(sigObj, score = "ideker")
        s <- (sumStat - bk$mu[2]) / bk$sigma[2]

        rn <- runif(iterations)
        rn2 <- runif(iterations)

        boundaries <- NULL

        # OPTIMISATION =========================================================

        for (i in seq_len(iterations)) {

            actNet <- sigObj@network[sigObj@network$active, ]$gene
            actNet <- as.character(actNet)
            adjSubgraph <- adjMatrix[actNet, ]

            if (length(dim(adjSubgraph))>0) {

                boundaries <- adjSubgraph[apply(adjSubgraph, 1, sum) > 0, ]
                boundaries <- colnames(
                    boundaries[, apply(boundaries, 2, sum) > 0]
                )
                # remove boundaries already active
                boundaries <- boundaries[!boundaries %in% actNet]

            } else {

                boundaries <- NULL

            }

            neighbours <- NULL
            probneighbour <- 0
            sampNeighbour <- FALSE
            lenActive <- length(actNet)

            if (lenActive > 2) {

                probneighbour <- 0.5

                if (probneighbour > rn2[i]) {

                    sampNeighbour <- TRUE

                    if (lenActive <= 10) {

                        CComp <- list(1, 2)

                        while (length(CComp) != 1) { #bottleneck

                            ge <- sample(actNet, 1)
                            test <- graph::subGraph(
                                as.character(actNet[actNet != ge]),
                                sigObj@connected_comp
                            )
                            CComp <- RBGL::connectedComp(test)

                            if (length(CComp) == 1) {
                                neighbours <- c(ge)
                            }

                        }

                    } else {

                        ge <- sample(actNet,1)
                        test <- graph::subGraph(
                            as.character(actNet[actNet != ge]),
                            sigObj@connected_comp
                        )
                        CComp <- RBGL::connectedComp(test)

                        if(max(sapply(CComp, length)) < 5) {
                            sampNeighbour <- FALSE
                        }

                        neighbours <- c(ge)

                    }
                }
            }

            if (lenActive > 1 & length(unique(c(neighbours, boundaries))) > 0) {

                # Sample a new gene to toggle its state
                if (!sampNeighbour) {
                    newG <- as.character(sample(unique(c(boundaries)), 1))
                } else {
                    newG <- neighbours
                }

                tog <- which(sigObj@network$gene %in% newG)
                sigObj@network[tog, ]$active <- !sigObj@network[tog, ]$active

                # Compute the subnet score
                if (!sampNeighbour) {

                    sumStat <- computeScore(sigObj, score = "ideker")

                    # Scale the subnet score
                    lenActive <- sum(sigObj@network$active)
                    s2 <- (sumStat - bk$mu[lenActive]) / bk$sigma[lenActive]

                } else {

                    if(length(CComp) == 1){

                        sc2 <- computeScore(sigObj, score = "ideker")
                        cK <- sum(sigObj@network$active)
                        s2 <- (sc2 - bk$mu[cK]) / bk$sigma[cK]

                    } else {

                        scs <- lapply(CComp, function(x) {
                            cK <- length(x)

                            if(cK >= 5) {
                                scn <- (1 / (sqrt(cK))) * sum(
                                    sigObj@network$score[
                                        as.character(sigObj@network$gene) %in% x
                                    ]
                                )
                                scn <- (scn - bk$mu[cK]) / bk$sigma[cK]
                            } else {
                                scn <- NA
                            }

                        }) #compute score for each CC

                        s2 <- max(unlist(scs), na.rm = TRUE)
                        GL <- CComp[[which(unlist(scs) == s2)[1]]]
                        nG <- as.character(sigObj@network$gene)
                        sigObj@network[!(nG %in% GL), ]$active <- FALSE

                    }
                }

                # Keep or not the toggled gene, acceptance probability
                if (s2 < s) {

                    # acceptance probability
                    prob <- exp((s2 - s) / sigObj@SA$temp[i])

                    if(i > iterations - iterations / 10) prob <- 0

                    if (prob > rn[i]) {

                        s <- s2
                        sigObj@SA$score[i] <- s
                        sigObj@SA$size[i] <- sum(sigObj@network$active)

                    } else {

                        cond <- which(sigObj@network$gene %in% newG)

                        sigObj@network[cond, ]$active  <-
                            !sigObj@network[cond,]$active

                        sigObj@SA$score[i] <- s
                        sigObj@SA$size[i] <- sum(sigObj@network$active)

                    }

                } else if (s2 >= s) {

                    s <- s2
                    sigObj@SA$score[i] <- s
                    sigObj@SA$size[i] <- sum(
                        sigObj@network$active
                    )

                }
            }
        }

        ### Return the results (subnetwork, size, score)
        sigObj@subnet_size <- sigObj@SA$size[iterations]

        sigObj@subnet_score <- sigObj@SA$score[iterations]

        sigObj@aggregate_score <-
            sum(sigObj@network$score[sigObj@network$active]) /
            sqrt(sigObj@subnet_size)

        sigObj@mean_score <-
            sum(sigObj@network$score[sigObj@network$active]) /
            (sigObj@subnet_size)

        sigObj@subnet_genes <- sigObj@network$gene[sigObj@network$active]

        invisible(sigObj)
    }
}
