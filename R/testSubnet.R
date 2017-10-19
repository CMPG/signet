#' Test the significance of high-scoring subnetworks found using simulated
#' annealing.
#'
#' @param sigObj A list of signet objects obtained using the
#' \verb{searchSubnet} function.
#' @param nullD Vector of null subnetwork scores generated using the
#' \code{nullDist} function.
#'
#' @return For each \verb{signet} object, a p-value is computed given the
#' provided emnpirical null distribution.
#'
#' @keywords null distribution, simulated annealing
#' @export
#' @examples
#' data(daub13)
#' \dontrun{
#' HSS <- searchSubnet(kegg_human[1:10], scores, iterations = 5000)
#' null <- nullDist(kegg_human, scores, n = 1000)
#' HSS <- testSubnet(HSS,null)
#' }

testSubnet<-function(sigObj, nullD) {

    if(class(sigObj) != "signetList") {

        stop("Input is not a list of signet objects")

    }

    if(class(nullD) != "numeric") {

        stop("Null distribution must be a numeric vector")

    }

    lapply(names(sigObj),function(x) {

        if(length(sigObj[[x]])>1) {

            sigObj[[x]]$p.value <<- mean(nullD > sigObj[[x]]$subnet_score,
                                         na.rm=TRUE)

        }

    })

    return(sigObj)

}
