#' @importFrom graphics axis layout mtext par plot
#' @importFrom stats runif sd
#' @importFrom utils installed.packages setTxtProgressBar txtProgressBar
#'
#' @keywords internal

computeScore <- function(signetObject, score = "ideker") {

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

    return(subnetStat)

}

#' @keywords internal
getAdjacencyMatrix <- function(pathway,
                               directed = FALSE,
                               selfLoops = FALSE) {

    requireNamespace("graph",quietly = TRUE)
    x <- graph::edges(pathway)
    GList <- names(x)
    adjMatrix <- matrix(0,
                        length(GList),
                        length(GList),
                        dimnames = list(GList,GList))

    if(!directed) {
        sapply(GList, function(y){
            adjMatrix[y,x[[y]]] <<-  1
            adjMatrix[x[[y]],y] <<-  1
        })
    } else {
        sapply(GList, function(y){
            adjMatrix[y,x[[y]]] <<-  1
        })
    }

    if(!selfLoops) diag(adjMatrix) <- 0

    rownames(adjMatrix) <- colnames(adjMatrix) <- GList

    return(adjMatrix)
}

#' @keywords internal
createSignetObject <- function(pathway, scores, iterations) {

    names(scores)[[1]] <- "gene";names(scores)[[2]] <- "score"
    scores <- scores[scores$gene %in% graph::nodes(pathway),]

    #to avoid duplicated nodes:
    X <- tapply(scores[!is.na(scores$score),]$score,
                scores[!is.na(scores$score),]$gene,max,na.rm=TRUE)

    X2 <- data.frame(gene=as.character(names(X)),score=X) ##num for entrez?

    scores <- rbind(X2,scores[is.na(scores$score) & !scores$gene %in% X2$gene,])

    connected_comp <- graph::subGraph(as.character(
        scores[!is.na(scores$score),]$gene),pathway)

    TA <- unlist(lapply(RBGL::connectedComp(connected_comp),length))
    maxi <- which(TA==max(TA))
    CC <- unlist(RBGL::connectedComp(connected_comp)[maxi])

    connected_comp <- graph::subGraph(CC[CC %in% as.character(
        scores[!is.na(scores$score),]$gene
    )],connected_comp)

    connected_comp <- graph::ugraph(connected_comp)

    nnodes <- dim(scores)[1]
    subnet_score <- NA
    subnet_size <- NA
    subnet_genes <- NA
    network <- data.frame(
        gene=scores$gene,
        score=scores$score,
        active=rep(FALSE,nnodes)
    )

    if(max(TA) < 5){

        object <- list(pathway=pathway,
                       connected_comp=connected_comp,
                       network=network,
                       subnet_score=subnet_score,
                       subnet_size=subnet_size,
                       subnet_genes=subnet_genes,
                       p.value=NA)
        class(object) <- "signet"
        return(object)

    } else {

        simulated_annealing <- data.frame(
            temperature=temperatureFunction(iterations),
            size_evolution=rep(NA,iterations),
            score_evolution=rep(NA,iterations)
        )

        object <- list(pathway=pathway,
                       connected_comp=connected_comp,
                       network=network,
                       subnet_score=subnet_score,
                       subnet_size=subnet_size,
                       subnet_genes=subnet_genes,
                       p.value=NA,
                       simulated_annealing=simulated_annealing)
        class(object) <- "signet"
        return(object)
    }
}

#' @export
print.signet <- function(x, ...) {

    cat("High-scoring subnetwork found with simulated annealing\n")

    cat(paste("Subnetwork score: ",
              round(x$subnet_score,digits=4),"\n",sep=""))

    cat(paste("Subnetwork size: ",
              x$subnet_size,"\n",sep=""))

    cat(paste("Genes in subnetwork: ",
              paste(x$subnet_genes,collapse=" "),"\n",sep=""))

}

#' @export
summary.signet <- function(object, ...) {

    cat("High-scoring subnetwork found with simulated annealing\n")

    cat(paste("Subnetwork score: ",
              round(object$subnet_score,digits=4),"\n",sep=""))

    cat(paste("Subnetwork size: ",
              object$subnet_size,"\n",sep=""))

    cat(paste("Genes in subnetwork: ",
              paste(object$subnet_genes,collapse=" "),"\n",sep=""))

}

#' @export
summary.signetList <- function(object, ...) {

    if(class(object) != "signetList") stop("Object must be a signetList")

    pathway <- names(object)

    net.size <- unlist(lapply(object, function(xl){
        if(length(xl)>1) return(dim(xl$network)[1])
        else return(NA)
    }))

    subnet.size <- unlist(lapply(object, function(xl){
        if(length(xl)>1) return(xl$subnet_size)
        else return(NA)
    }))

    subnet.score <- unlist(lapply(object, function(xl){
        if(length(xl)>1) return(xl$subnet_score)
        else return(NA)
    }))

    p.value <- unlist(lapply(object, function(xl){
        if(length(xl)>1) return(xl$p.value)
        else return(NA)
    }))

    subnet.genes <- unlist(lapply(object,function(xl){
        if(length(xl)>1) return(paste(xl$subnet_genes,collapse=" "))
        else return(NA)
    }))

    signet_table <- data.frame(pathway,
                               net.size,
                               subnet.size,
                               subnet.score,
                               p.value,
                               subnet.genes)

    rownames(signet_table) <- NULL
    return(signet_table)
}

#' @export
plot.signet <- function(x, ...) {

    if(class(x) != "signet") {
        stop("Input must be a signet object.")
    }

    if(is.na(x$subnet_score)) {
        stop("Pathway is too small/disconnected to be tested.")
    }

    requireNamespace("graphics", quietly = TRUE)

    m <- rbind(c(0,1,1,1,1,0),c(2,2,2,3,3,3))
    layout(m)

    glist <- x$network[x$network$active,]$gene

    subs <- as.character(glist)
    subg <- graph::subGraph(subs,x$connected_comp)

    col <- c(rep("red",length(subs)))
    nAttrs <- list()
    nAttrs$fillcolor <- col
    nAttrs$height <- rep("0.6", length(graph::nodes(x$connected_comp)))
    nAttrs$width <- rep("0.6", length(graph::nodes(x$connected_comp)))
    nAttrs$color <- rep("darkgrey", length(graph::nodes(x$connected_comp)))

    names(nAttrs$color) <- names(nAttrs$width) <- names(nAttrs$height) <-
        graph::nodes(x$connected_comp)

    names(nAttrs$fillcolor) <- c(subs)

    eAttrs <- list()
    eAttrs$color <- rep("grey",length(graph::edgeNames(x$connected_comp)))
    names(eAttrs$color) <- graph::edgeNames(x$connected_comp)

    graph::plot(x$connected_comp, y = "neato", nodeAttrs = nAttrs,
                edgeAttrs = eAttrs)

    x.ax <- 1:length(x$simulated_annealing$temperature)
    y <- x$simulated_annealing$size_evolution
    z <- x$simulated_annealing$score_evolution
    par(mar = c(5, 5, 5, 5))  # Leave space for z axis

    plot(x.ax, z, type = "l",lwd = 1,cex = 0.2,pch = 16,
         col = "firebrick",ylab = "Subnetwork score",
         xlab = "Iterations") # first plot

    par(new = TRUE)
    plot(x.ax, x$simulated_annealing$temperature, type = "l", lwd=1,
         xlab = "", ylab = "", axes=FALSE, lty=2,
         col = "grey")
    axis(side = 4, at = pretty(range(x$simulated_annealing$temperature)))
    mtext("Temperature", side = 4, line = 3, cex=0.7)

    plot(x.ax, y, type = "l", cex = 0.5,
         xlab = "Iterations", ylab = "Subnetwork size",
         pch=16,lwd=1,col="dodgerblue")

    par(new=TRUE)
    plot(x.ax, x$simulated_annealing$temperature, type="l",lwd=1,
         xlab="", ylab="",axes=FALSE,lty=2,
         col="grey")
    axis(side=4, at = pretty(range(x$simulated_annealing$temperature)))
    mtext("Temperature", side=4, line=3, cex=0.7)

    par(mfrow = c(1,1))

}


#' @keywords internal
adjacencyMatrixToList <- function(adjMatrix) {

    adjList <- apply(adjMatrix, 1, function(x) return(names(x[x==1])))
    return(adjList)

}

#' @keywords internal
lapply_pb <- function(X, FUN, ...) {

    env <- environment()
    pb_Total <- length(X)
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

    # wrapper around FUN
    wrapper <- function(...) {
        curVal <- get("counter", envir = env)
        assign("counter", curVal +1 ,envir=env)
        setTxtProgressBar(get("pb", envir=env), curVal +1)
        FUN(...)
    }

    res <- lapply(X, wrapper, ...)
    close(pb)
    res

}

#' @keywords internal
temperatureFunction <- function(iterations, threshold=0.02) {

    alpha <- exp(log(threshold)/iterations)
    t <- 1*alpha^(1:iterations)
    return(t)

}
