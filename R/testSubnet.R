#' Test the significance of high-scoring subnetworks
#'
#' @param sigObj A list of signet objects.
#' @param nullD Vector of null highscores generated using the
#' \code{nullDist} function
#'
#' @keywords null distribution, simulated annealing
#' @export


testSubnet<-function(sigObj,nullD) {

  if(class(sigObj) != "signetList") {

    stop("Input is not a list of signet objects")

  }

  if(class(nullD) != "numeric") {

    stop("Null distribution must be a numeric vector")

  }

  lapply(names(sigObj),function(x) {

    if(length(sigObj[[x]])>1) {

      sigObj[[x]]$p.value <<- mean(nullD > sigObj[[x]]$subnet_score, na.rm=TRUE)

    }

  })

  return(sigObj)

}
