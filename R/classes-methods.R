#' @importFrom graph graphNEL
#' @importFrom methods new show initialize
NULL

#' An S4 class to represent a pathway and the results of the associated
#' simulated annealing run.
#'
#'
#' @slot connected_comp A graphNEL object (biological pathway)
#' @slot network A data frame (gene IDs and scores)
#' @slot SA A data frame (information on the simulated annealing run)
#' @slot subnet_score A numeric value (subnetwork score)
#' @slot aggregate_score A numeric value (aggregate subnetwork score)
#' @slot mean_score A numeric value (average gene score in the pathway)
#' @slot subnet_size An integer value (subnetwork size)
#' @slot subnet_genes A factor (subnetwork genes)
#' @slot p.value A numeric value (empirical p-value)
#'
#' @name Signet-class
#' @rdname Signet-class
#' @export
.Signet <- setClass(
    "Signet",
    slots = c(
        connected_comp = "graphNEL",
        network = "data.frame",
        SA = "data.frame",
        subnet_score = "numeric",
        aggregate_score = "numeric",
        mean_score = "numeric",
        subnet_size = "integer",
        subnet_genes = "factor",
        p.value = "numeric"
    )
)

#' @describeIn Signet Print the summary a Signet object
#' @param object A signet object.
#' @export
setMethod("show", "Signet", function(object) {
    cat(paste0(
        "Subnetwork score: ", round(object@subnet_score, digits = 4), "\n")
    )
    cat(paste0("Subnetwork size: ", object@subnet_size, "\n"))
    cat(paste0(
        "Genes in subnetwork: ", paste(
            object@subnet_genes, collapse = " "
        ), "\n"
    ))
})

#' @describeIn Signet Print the summary of a Signet object
#' @export
setMethod("summary", "Signet", function(object) {
    show(object)
})

#' @describeIn Signet Plot a Signet object
#' @param x A signet object.
#' @param y Omitted when plotting a Signet object.
#' @param ... Other graphical parameters.
#' @return A plot of the simulated annealing run.
#' @export
setMethod("plot", c("Signet", "missing"), function(x, y, ...) {

    if(is.na(x@subnet_score)) {
        stop("Pathway is too small or disconnected to be tested.")
    }

    requireNamespace("graphics", quietly = TRUE)

    par(mfrow = c(1, 2))

    x.ax <- seq_len(length(x@SA$temp))
    y <- x@SA$size
    z <- x@SA$score
    par(mar = rep(5, 4))  # Leave space for z axis

    plot(
        x.ax, z, type = "l", lwd = 1, cex = 0.2, pch = 16,
        col = "firebrick", ylab = "Subnetwork score", xlab = "Iterations"
    )

    par(new = TRUE)
    plot(
        x.ax, x@SA$temp,
        type = "l", lwd = 1, xlab = "", ylab = "", axes = FALSE, lty = 2,
        col = "grey"
    )
    axis(side = 4, at = pretty(range(x@SA$temp)))
    mtext("Temperature", side = 4, line = 3, cex = 1)

    plot(
        x.ax, y, type = "l", cex = 0.5,
        xlab = "Iterations", ylab = "Subnetwork size",
        pch = 16,lwd = 1,col = "dodgerblue"
    )

    par(new = TRUE)
    plot(
        x.ax, x@SA$temp,
        type = "l", lwd = 1, xlab = "", ylab = "", axes = FALSE, lty = 2,
        col = "grey"
    )
    axis(side = 4, at = pretty(range(x@SA$temp)))
    mtext("Temperature", side = 4, line = 3, cex = 1)

    par(mfrow = c(1,1))

})

#' @describeIn Signet Initialize a Signet object
#' @param .Object Object to initialize.
#' @param pathway Biological pathway (graphNEL object).
#' @param scores Gene scores list.
#' @param iterations Number of simulated annealing iterations.
#' @return A signet object.
#' @export
setMethod("initialize", "Signet",
        function(.Object, pathway, scores, iterations) {

    names(scores)[[1]] <- "gene"
    names(scores)[[2]] <- "score"
    scores <- scores[scores$gene %in% graph::nodes(pathway),]

    #to avoid duplicated nodes:
    X <- tapply(
        scores[!is.na(scores$score),]$score,
        scores[!is.na(scores$score),]$gene,
        max,
        na.rm = TRUE
    )

    X2 <- data.frame(gene = as.character(names(X)), score = X)

    scores <- rbind(
        X2,
        scores[is.na(scores$score) & !scores$gene %in% X2$gene,]
    )

    geList <- as.character(scores[!is.na(scores$score),]$gene)

    connected_comp <- graph::subGraph(geList, pathway)

    TA <- unlist(lapply(RBGL::connectedComp(connected_comp), length))
    maxi <- which(TA == max(TA))
    CC <- unlist(RBGL::connectedComp(connected_comp)[maxi])

    connected_comp <- graph::subGraph(CC[CC %in% geList], connected_comp)

    .Object@connected_comp <- graph::ugraph(connected_comp)

    nnodes <- dim(scores)[1]

    .Object@subnet_score <- as.numeric(NA)
    .Object@subnet_size <- as.integer(NA)
    .Object@subnet_genes <- as.factor(NA)
    .Object@network <- data.frame(
        gene = scores$gene,
        score = scores$score,
        active = rep(FALSE, nnodes)
    )

    .Object@SA <- data.frame(
        temp = temperatureFunction(iterations),
        size = rep(NA, iterations),
        score = rep(NA, iterations)
    )

    .Object@p.value = as.numeric(NA)

    .Object
})



#' An S4 class to represent a list of "Signet" objects.
#'
#' @return Results of the simulated annealing run for a lisst of pathways.
#' @slot results A list of Signet objects.
#' @name SignetList-class
#' @rdname SignetList-class
#' @export
.SignetList <- setClass(
    "SignetList",
    slots = c(
        results = "list"
    )
)

#' @describeIn SignetList Initialize a SignetList
#' @param .Object A SignetList object.
#' @param list A list of Signet objects.
#' @return A SignetList object.
#' @export
setMethod("initialize", "SignetList", function(.Object, list) {

    .Object@results <- list
    .Object

})

#' @describeIn SignetList Access the ith element (signet object) of the
#' SignetList
#' @param x A SignetList object.
#' @param i Index specifying elements to extract or replace.
#' @aliases [[
#' @export
setMethod("[[", "SignetList", function(x, i) {

    x@results[[i]]

})

#' @describeIn SignetList Access the ith element (signet object) of the
#' SignetList
#' @param j Unused for SignetList objects.
#' @param drop Unused for SignetList objects.
#' @param ... Unused for SignetList objects.
#' @export
setMethod("[", c("SignetList", "ANY", "missing"),
        function(x, i, j, ..., drop = TRUE)  {

    x@results[i]

})

#' @describeIn SignetList Summarize the SignetList.
#' @param object A SignetList object.
#' @return A data frame containing
#' summary statistics for each element (network and subnetwork sizes,
#' subnetwork score, p-value, significant genes list)
#' @export
setMethod("summary", "SignetList", function(object, ...) {

    re <- data.frame(pathway = names(object@results))

    xu <- sapply(as.character(re$pathway), function(xl) {
        pa <- object@results[[xl]]
        return(c(
            net.size = dim(pa@network)[1],
            subnet.size = pa@subnet_size,
            subnet.score = pa@subnet_score,
            p.val = pa@p.value,
            subnet.genes = paste(pa@subnet_genes, collapse = " ")
        ))
    })

    re <- cbind(re, t(xu))
    rownames(re) <- NULL

    re

})
