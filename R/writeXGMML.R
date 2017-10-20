#' Write Cytoscape input file
#'
#' This function allows to write an XGMML file to represent the results in
#' Cytoscape.
#'
#' @param results A signet or signetList object.
#' @param filename The desired file name. Default is "signet_output.xgmml".
#' @param threshold Significance threshold (default: 0.01). If a signetList is
#' provided, all subnetworks with a p-value below this threshold will be merged
#' and represented.
#'
#' @keywords subnetwork, network, plot, Cytoscape, visualization
#'
#' @return Writes an XGMML file in the working directory. If a single pathway
#' (signet object) is provided, the whole pathway is represented and nodes
#' belonging to the highest-scoring subnetwork (HSS) are highlighted in red. If
#' a list of pathways (signetList) is provided, all subnetworks with a p-value
#' below a given threshold (default: 0.01) are merged and represented. Note that
#' in this case, only the nodes belonging to HSS are kept for representation.
#'
#' @export
#' @examples
#' data(daub2013) #pathways and gene scores from Daub et al. (2013).
#' \dontrun{
#' #run the search in all the pathways with 5000 iterations (default)
#' example <- searchSubnet(kegg_human[[1]], scores)
#' writeXGMML(example)
#' }

writeXGMML <- function(results,
                       filename = "signet_output.xgmml",
                       threshold = 0.01) {

    if (sum(installed.packages()[, 1]=="igraph")==0) {
        stop("Package igraph is required.")
    } else {
        requireNamespace("igraph", quietly=TRUE)
    }


    if (class(results) != "signet" & class(results) != "signetList") {
        stop("results must be a signet or signetList object.")
    }

    if(substr(filename, nchar(filename)-5, nchar(filename)) != ".xgmml") {
        filename <- paste0(filename, ".xgmml")
    }

    if (class(results) == "signet") {

        el <- c()
        for (i in 1:length(graph::edges(results$connected_comp))) {
            el <- c(el, paste(
                names(graph::edges(results$connected_comp)[i]),
                graph::edges(results$connected_comp)[[i]],
                sep = "-"
            ))
        }
        xx <- do.call("rbind", strsplit(el, "-"))

        for (i in 1:nrow(xx)) {
            xx[i,] = sort(xx[i,])
        }

        xx = xx[!duplicated(xx), ]

        source <- xx[, 1]
        target <- xx[, 2]

        el <- apply(xx, 1, function(x)
            (paste(x, collapse = "-")))

        nodelab <- results$network$gene
        nodescores <- results$network$score
        nodeact <- results$network$active

        nsize <-
            35 + ((80 - 35) *
                      (results$network$score - min(results$network$score))) /
            (-min(results$network$score) + max(results$network$score))

        ncol <- rep(NA, length(nsize))
        ncol[results$network$active] <- "#CC0000"
        ncol[!results$network$active] <- "#FFFFFF"

        gr <- igraph:::graph_from_graphnel(results$connected_comp)
        ncoords <- igraph::layout.fruchterman.reingold(gr) * 50

    } else {

        selec <- unlist(lapply(results, function(x) {
            if (is.na(x$p.value)) return(FALSE)
            if (length(x) > 1)
                return(x$p.value < threshold)
            else
                return(FALSE)
        }))
        results <- results[selec]

        if(length(results) == 0) {
            stop("No subnetwork has a p-value lower than the specified
                 threshold.")
        }

        subG <-
            graph::subGraph(as.character(results[[1]]$subnet_genes),
                            results[[1]]$connected_comp)
        GRAPH <- igraph:::graph_from_graphnel(subG)

        # merge graphs
        for (i in 2:length(results)) {
            subG <- graph::subGraph(as.character(results[[i]]$subnet_genes),
                                    results[[i]]$connected_comp)
            subG <- igraph:::graph_from_graphnel(subG)
            GRAPH <- igraph::graph.union(GRAPH, subG)
        }

        newObj <- list()
        newObj$connected_comp <- igraph:::as_graphnel(GRAPH)

        newObj$network <- do.call("rbind", lapply(results, function(x) {
            x$network[x$network$active,]
        }))

        nodelab <- newObj$network$gene
        nodescores <- newObj$network$score
        nodeact <- newObj$network$active


        el <- c()
        for (i in 1:length(graph::edges(newObj$connected_comp))) {
            el <- c(el, paste(
                names(graph::edges(newObj$connected_comp)[i]),
                graph::edges(newObj$connected_comp)[[i]],
                sep = "-"
            ))
        }
        xx <- do.call("rbind", strsplit(el, "-"))

        for (i in 1:nrow(xx)) {
            xx[i,] = sort(xx[i,])
        }

        xx = xx[!duplicated(xx), ]

        source <- xx[, 1]
        target <- xx[, 2]

        el <- apply(xx, 1, function(x)
            (paste(x, collapse = "-")))
        nsize <- 35 + ((80 - 35) * (nodescores - min(nodescores))) /
            (-min(nodescores) + max(nodescores))

        ncol <- rep("#CC0000", length(nsize))
        ncoords <- igraph::layout.fruchterman.reingold(GRAPH) * 50

    }

    cat(
        paste(
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
            <graph label="signet results"
            xmlns:dc="http://purl.org/dc/elements/1.1/"
            xmlns:xlink="http://www.w3.org/1999/xlink"
            xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
            xmlns:cy="http://www.cytoscape.org"
            xmlns="http://www.cs.rpi.edu/XGMML"
            directed="0">'
        ),

        paste(
            '\n\n<node label="',
            nodelab,
            '" id="',
            nodelab,
            '">
            <att name="size" type="real" value="',
            nodescores,
            '"/>
            <att name="subnet" type="boolean" value="',
            nodeact,
            '"/>
            <graphics z="0.0" outline="#000000" x="',
            ncoords[, 1],
            '"  type="ELLIPSE" w="',
            nsize,
            '" width="2.0" fill="',
            ncol,
            '" h="',
            nsize,
            '" y="',
            ncoords[, 2],
            '" >
            <att name="NODE_TOOLTIP" value="" type="string"/>
            <att name="NODE_LABEL_FONT_FACE"
                value="SansSerif.plain,plain,10" type="string"/>
            <att name="NODE_VISIBLE" value="true" type="string"/>
            <att name="NODE_BORDER_STROKE" value="SOLID" type="string"/>
            <att name="NODE_CUSTOMGRAPHICS_SIZE_8"
                value="41.53651611954727" type="string"/>
            <att name="NODE_LABEL" value="',
            nodelab,
            '" type="string"/>
            <att name="NODE_LABEL_COLOR" value="#FFFFFF" type="string"/>
            <att name="NODE_DEPTH" value="0.0" type="string"/>
            <att name="NODE_NESTED_NETWORK_IMAGE_VISIBLE"
                value="true" type="string"/>
            <att name="NODE_LABEL_FONT_SIZE" value="10" type="string"/>
            <att name="NODE_TRANSPARENCY" value="255" type="string"/>
            <att name="NODE_BORDER_TRANSPARENCY" value="255" type="string"/>
            <att name="NODE_LABEL_TRANSPARENCY" value="255" type="string"/>
            </graphics>
            </node>',
            sep = ''
        ),

        paste(
            '\n\n<edge label="',
            el,
            '" source="',
            source,
            '" target="',
            target,
            '">
            <att name="weight" type="integer" value="1"/>
            <graphics width="1.5" fill="#000000">
            <att name="EDGE_TARGET_ARROW_SHAPE" value="NONE" type="string"/>
            </graphics>
            </edge>',
            sep = ''
        ),

        paste('\n</graph>'),
        file = filename
    )
}
