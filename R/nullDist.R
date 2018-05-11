#' Null distribution
#'
#' Generate the high-scores null distribution to compute empirical p-values
#' for each biological pathway.
#'
#' @param pathways A list of graphNEL objects.
#' @param scores A data frame in which the first column corresponds to the gene
#' ID and the second column contains the gene scores.
#' @param n Number of null high-scores to compute (default = 1000).
#' @param background Optional. Background distribution computed using the
#' \verb{backgroundDist} function.
#'
#' @return A vector of subnetworks scores obtained under the null hypothesis.
#' Must be used as input of the \verb{testSubnet} function.
#'
#' @keywords simulated annealing, null distribution
#' @export
#' @examples
#' # Get KEGG pathways from the package graphite:
#' # library(graphite)
#' # kegg <- pathways("hsapiens", "kegg")
#' # kegg_human <- lapply(kegg, pathwayGraph)
#'
#' data(daub13) # load the gene scores
#'
#' # generate the null distribution (here, only 5 values, but
#' # at least 1000 are advised)
#' null <- nullDist(kegg_human, scores, n = 5)

nullDist <- function(pathways,
                        scores,
                        n = 1000,
                        background) {

    requireNamespace("graph", quietly = TRUE)
    colnames(scores) <- c("gene", "score")

    if (missing(background)) {
        sizeP <- 200
        background <- data.frame(
            mu = sqrt(seq_len(sizeP)) * mean(scores$score, na.rm = TRUE),
            sigma = rep(sd(scores$score, na.rm = TRUE), sizeP)
        )
    }

    genesDB <- unique(unlist(sapply(pathways, graph::nodes)))
    scores <- scores[scores$gene %in% genesDB, ]

    message("  Generating null distribution...")

    bk <- lapply_pb(seq_len(n), function(x) {

        newscores <- data.frame(
            gene = scores$gene,
            score = sample(scores$score)
        )

        cond <- TRUE

        while(cond) {

            grap <- sample(seq_len(length(pathways)), 1)

            HSS <- try(searchSubnet(
                pathways[[grap]],
                scores = newscores,
                background,
                iterations = 1000
            ), silent=TRUE)

            if(class(HSS)!="try-error") cond <- FALSE

        }

        return(HSS@subnet_score)

    })

    return(unlist(bk))
}
