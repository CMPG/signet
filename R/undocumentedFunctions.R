#' @importFrom graphics axis layout mtext par plot
#' @importFrom stats runif sd
#' @importFrom utils installed.packages setTxtProgressBar txtProgressBar

#' @keywords internal
computeScore <- function(signetObject, score = "ideker") {

    activeNet <- signetObject@network[signetObject@network$active, ]

    if (score ==  "mean") subnetStat <- mean(activeNet$score)
    if (score ==  "sum") subnetStat <- sum(activeNet$score)
    if (score ==  "ideker") {
        k <- length(activeNet$score)
        subnetStat <- (1 / sqrt(k)) * sum(activeNet$score)
    }

    return(subnetStat)

}

#' @keywords internal
getAdjacencyMatrix <- function(pathway, directed = FALSE, selfLoops = FALSE) {

    requireNamespace("graph", quietly = TRUE)
    x <- graph::edges(pathway)
    GList <- names(x)

    adjMatrix <- matrix(
        0,
        length(GList),
        length(GList),
        dimnames = list(GList, GList)
    )

    e <- environment()
    if(!directed) {
        sapply(GList, function(y) {
            tmp <- get("adjMatrix", envir = e)
            tmp[y, x[[y]]] <- 1
            tmp[x[[y]], y] <- 1
            assign(x = "adjMatrix", value = tmp, envir = e)
        })
    } else {
        sapply(GList, function(y) {
            tmp <- get("adjMatrix", envir = e)
            tmp[y, x[[y]]] <- 1
            assign(x = "adjMatrix", value = tmp, envir = e)
        })
    }

    if(!selfLoops) diag(adjMatrix) <- 0

    rownames(adjMatrix) <- colnames(adjMatrix) <- GList

    return(adjMatrix)
}

#' @keywords internal
adjacencyMatrixToList <- function(adjMatrix) {

    adjList <- apply(adjMatrix, 1, function(x) {
        names(x[x ==  1])
    })

    return(adjList)

}

#' @keywords internal
lapply_pb <- function(X, FUN, ...) {

    env <- environment()
    pb_Total <- length(X)
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

    wrapper <- function(...) {
        setTxtProgressBar(pb, counter)
        counter <<- counter + 1L
        FUN(...)
    }

    res <- lapply(X, wrapper, ...)
    close(pb)
    res

}

#' @keywords internal
temperatureFunction <- function(iterations, threshold = 0.02) {

    alpha <- exp(log(threshold) / iterations)
    t <- 1 * alpha ^ seq_len(iterations)
    return(t)

}
