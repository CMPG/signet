#' Temperature function
#'
#' Generate the temperature function used in simulated annealing.
#' @param iterations A graph object.
#' @param param a data frame with gene list and associated scores
#' @param burnin
#' @keywords subnetwork, simulated annealing
#' @export
#' @examples
#' #get data
#' temperatureFunction()

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
