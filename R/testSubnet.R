#' Test the significance of a subnetwork
#'
#' @keywords null distribution, simulated annealing
#' @export


testSubnet<-function(sigObj,nullD) {

  if(class(sigObj) != "signetList") stop("Input is not a list of signet objects")
  if(class(tt) != "numeric") stop("Null distribution must be a numeric vector")

  lapply(names(sigObj),function(x) {

    if(length(sigObj[[x]])>1) {
      sigObj[[x]]$p.value <<- mean(tt > sigObj[[x]]$subnet_score, na.rm=TRUE)
    }

  })

  return(sigObj)

}
