#' Generate the temperature function
#'
#' A function to generate the temperature function used in simulated annealing.
#' After a burnin period, the temperature decreases geometrically.
#'
#' @param iterations Number of iterations.
#' @param param Geometric parameter
#' @param burnin Burnin
#' @keywords subnetwork, simulated annealing
#' @export
#' @examples
#' t<-temperatureFunction(iterations=1500,param=0.995,burnin=100)
#' plot(t)

temperatureFunction<-function(iterations,threshold=0.01)
{
  alpha <- exp(log(threshold)/iterations)
  t<-array(0,iterations)
  t[1]<-1
  for(i in 2:iterations)
  {
    t[i]<-alpha*t[i-1]
  }
  return(t)
}
