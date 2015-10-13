#' Test the significance of a subnetwork
#'
#' A function to generate the temperature function used in simulated annealing.
#' After a burnin period, the temperature decreases geometrically.
#'
#' @param outputSA Output of the algorithm.
#' @param pathwaysList List of pathways to generate null distribution.
#'
#' @keywords subnetwork, simulated annealing
#' @export
#' @examples
#' #t<-temperatureFunction(iterations=1500,param=0.995,burnin=100)
#' #plot(t)

testSubnet<-function(outputSA,pathwaysList,scores)
{

  nul<-nullDistribution(pathwaysList,scores,kmin=outputSA$size,kmax=outputSA$size,
                   iterations=10000)

  normScore<-(outputSA$score-nul$mu)/nul$sigma
  pval<-pnorm(normScore)

  cat("\n  p-value =",pval,"\n")

  outputSA$pvalue<-pval
  invisible(outputSA)
}
