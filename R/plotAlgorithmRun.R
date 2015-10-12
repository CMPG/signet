#' Generate an animation representing the subnetwork research
#'
#' This function will create an HTML file to visualize a run.
#' @param pathway A graph object.
#' @param scores a data frame with gene list and associated scores
#' @param nullDist A data frame with three columns : k, mu, sigma. Can be obtained thanks to nulldistribution() function
#' @param iterations The number of iterations for simulated annealings
#' @param temperature The temperature function parameter.
#' @param kmin The minimal size of a subnetwork.
#' @param directed If TRUE, considers the edges direction, i.e. cannot go back.
#' @param verbose If TRUE, displays text in the R console.
#' @param animPlot Number of iterations represented in the animation.
#' @keywords subnetwork, simulated annealing, heuristics, search algorithm
#' @return A HTML file and three folder in the animation subfolder in your working directory.
#' @export
#' @examples
#' require(signet)
#' data(keggPathways)
#' data(zScores)
#'
#' #plotAlgorithmRun(keggPathways[[1]],zScores,animPlot=10)

plotAlgorithmRun<-function(pathway,
                           scores,
                           nullDist,
                           iterations = 1000,
                           temperature = 0.990,
                           kmin = 5,
                           directed = FALSE,
                           verbose = TRUE,
                           animPlot = 100)
{
  requireNamespace("animation",quietly=TRUE)
  animation::saveHTML({
  searchSubnet(pathway,
               scores,
               nullDist,
               iterations,
               temperature,
               kmin,
               directed,
               verbose,
               animPlot)
  }, img.name = "algorithm_plot", imgdir = "unif_dir", htmlfile = "algorithm.html",
  autobrowse = FALSE, title = "Demo of simulated annealing",
  description = "Description",ani.height = 800, ani.width = 800)

  cat("\nDone !\n")
}
