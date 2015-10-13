#' Test the significance of a subnetwork
#'
#' A function to generate the temperature function used in simulated annealing.
#' After a burnin period, the temperature decreases geometrically.
#'
#' @param outputSA Output of the algorithm.
#' @param pathwaysList List of pathways to generate null distribution.
#' @param scores Again, the scores.
#'
#' @keywords subnetwork, simulated annealing
#' @export
#' @examples
#' #t<-temperatureFunction(iterations=1500,param=0.995,burnin=100)
#' #plot(t)

testSubnet<-function(outputSA,pathwaysList,scores)
{

  klist<-unlist(lapply(outputSA,function(x) return(x$size)))

  nul<-nullDistribution(pathwaysList,scores,kmin=min(klist),kmax=max(klist),
                   iterations=100)

  pvals<-lapply(outputSA,
  function(x)
  {
    normScore<-(x$score-nul[which(nul$k==x$size),]$mu)/nul[which(nul$k==x$size),]$sigma
    pval<-1-pnorm(normScore)
    outputSA$pvalue<-pval
    print(pval)
  })

  invisible(outputSA)
}
