#' Null distribution
#'
#' Generate the temperature function used in simulated annealing.
#' @param iterations A graph object.
#' @param param a data frame with gene list and associated scores
#' @param burnin A burnin
#' @keywords subnetwork, simulated annealing
#' @export
#' @examples
#' #get data
#' temperatureFunction()

nullDistribution<-function(pathwaysList,
                           scores,
                           kmin,
                           kmax,
                           iterations)
{
  require(graph)
  if(missing(pathwaysList))
  {
    pathwaysList<-NULL
    gList<-scores
  }
  cat("  Generating null distribution...\n")

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
          glis<-nodes(path)
        }
        gList<-scores[scores$gene %in% glis,]
      }
      ba<-rbind(ba,c(mean(gList[sample(length(gList$gene),k),]$score),na.rm=TRUE))
    }
    nullD<-rbind(nullD,c(k,mean(ba,na.rm=TRUE),sd(ba,na.rm=TRUE)))
  }
  cat("\n  Done !\n\n")
  nullD<-as.data.frame(nullD)
  colnames(nullD)<-c("k","mu","sigma")
  # print(nullD)
  return(nullD)
}
