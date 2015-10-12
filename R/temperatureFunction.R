#' Temperature function
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

temperatureFunction<-function(iterations,param,burnin)
{
  if(missing(burnin)) burnin <- iterations/20

  t<-array(0,iterations)
  t[1:burnin]<-1
  for(i in burnin:iterations)
  {
    t[i]<-param*t[i-1]
  }
  # t[t<1e-3]<-0
  return(t)
}
