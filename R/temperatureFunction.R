#' Generate the temperature function
#'
#' A function to generate the temperature function used in simulated annealing.
#' The temperature starts from 1 and decreases geometrically for a given number
#' of iterations.
#'
#' @param iterations Number of iterations.
#' @param threshold Temperature desired at the last iteration.
#' @keywords subnetwork, simulated annealing, temperature
#' @export
#' @examples
#' t <- temperatureFunction(iterations = 5000)
#' plot(t, ylab = "Temperature", xlab = "Iterations")

temperatureFunction<-function(iterations,threshold=0.02)
{
  alpha <- exp(log(threshold)/iterations)
  t<-1*alpha^(1:iterations)
  return(t)
}
