#' Test the significance of high-scoring subnetworks found using simulated
#' annealing.
#'
#' @param sigObj A list of signet objects obtained using the
#' \verb{searchSubnet} function.
#' @param null Vector of null subnetwork scores generated using the
#' \code{nullDist} function.
#'
#' @return For each \verb{signet} object, a p-value is computed given the
#' provided emnpirical null distribution.
#'
#' @keywords null distribution, simulated annealing
#' @export
#' @examples
#'
#' # Get KEGG pathways from the package graphite:
#' library(graphite)
#' kegg <- pathways("hsapiens", "kegg")
#' kegg_human <- lapply(kegg[1:5], pathwayGraph)
#'
#' data(daub13) # load the gene scores from Daub et al. (2013)
#'
#' #run the search in all the pathways with 2500 iterations (default)
#' example <- searchSubnet(kegg_human, scores)
#'
#' # generate the null distribution (here, only 10 values, but
#' # at least 1000 are advised)
#' null <- nullDist(kegg_human, scores, n = 10)
#' example <- testSubnet(example, null) #now, 'example' includes p-values
#' summary(example)

testSubnet <- function(sigObj, null) {

    if(class(sigObj) != "SignetList") {
        stop("Input is not a list of signet objects")
    }

    if(class(null) != "numeric") {
        stop("Null distribution must be a numeric vector")
    }

    sigObj@results <- lapply(sigObj@results, function(x) {
        x@p.value <- mean(null > x@subnet_score, na.rm = TRUE)
        x
    })

    return(sigObj)

}
