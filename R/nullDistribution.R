#' Null distribution
#'
#' Generate the temperature function used in simulated annealing.
#' @param pathwaysList A graph object, or a liste of graphNEL objects.
#' @param scores A data frame with gene list and associated scores
#' @param kmin kmin
#' @param kmax kmax
#' @param iterations Number of iterations.
#' @param distribution If TRUE, returns the whole simulated data. Else, returns the mean and sd.
#'
#' @keywords subnetwork, simulated annealing
#' @export
#' @examples
#' require(signet)
#' #get data
#' #nullDistribution()

nullDistribution<-function(pathwaysList,
                           scores,
                           kmin,
                           kmax,
                           iterations,
                           distribution=FALSE)
{
  requireNamespace("graph",quietly=TRUE)
  if(missing(pathwaysList))
  {
    pathwaysList<-NULL
    gList<-scores
  }
 colnames(scores)<-c("gene","score")

  cat("  Generating null distribution...\n")
  distrib<-NULL
  nullD<-NULL
  for(k in kmin:kmax)
  {
    # flush.console()
    cat("\r  ... for k =",k)

    ba<-NULL
    for(i in 1:iterations)
    {
      if(length(pathwaysList)>0)
      {

        glis<-NULL
        while(length(glis)<k+5)
        {
          path<-pathwaysList[[sample(length(pathwaysList),1)]]
          glis<-graph::nodes(path)
        }
        gList<-scores[scores$gene %in% glis,]
      }
      ba<-rbind(ba,c(mean(gList[sample(length(gList$gene),k,replace=TRUE),]$score,na.rm=TRUE)))
    }

    nullD<-rbind(nullD,c(k,mean(ba,na.rm=TRUE),sd(ba,na.rm=TRUE)))
    distrib<-c(distrib,ba)
  }
  cat("\n  Done !\n\n")
  nullD<-as.data.frame(nullD)
  colnames(nullD)<-c("k","mu","sigma")

  if(distribution) nullD<-list(parameters=nullD,distributions=distrib)
  # print(nullD)
  return(nullD)
}
