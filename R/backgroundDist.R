#' Background distribution
#'
#' Generate the background distribution necessary to compute subnetwork scores
#' during a simulated annealing run.
#'
#' @param pathwaysList A graphNEL object, or a list of graphNEL objects.
#' @param scores A data frame in which the first column if the gene list
#' and the second column contains the gene scores.
#' @param kmin Minimal value of k for which a background distribution is
#' generated. Default value is 1.
#' @param kmax Maximal value of k for which a background distribution is
#' generated. Default value is the size of the biggest graph in the
#' provided list.
#' @param iterations Number of iterations to make the bachground distribution.
#' @param distribution If TRUE, returns the whole simulated data.
#' Else, returns the mean and standard deviation of the distribution for each k.
#'
#' @keywords subnetwork, simulated annealing
#' @export
#' @examples
#' require(signet)
#' #get data
#' #backgroundDist()

backgroundDist<-function(pathwaysList,
                         scores,
                         kmin = 1,
                         kmax,
                         iterations = 1000,
                         subnetScore="sum",
                         distribution = FALSE)
{
  requireNamespace("graph",quietly=TRUE)
  if(missing(pathwaysList))
  {
    pathwaysList<-NULL
    gList<-scores
  }
 colnames(scores)<-c("gene","score")

 if(missing(kmax)) # kmax not specified: max k in the graph list
 {
   kmax <- max(unlist(lapply(pathwaysList,function(x)length(nodes(x)))))
 }

  cat("  Generating null distribution...\n")
  distrib<-NULL
  nullD<-NULL
  for(k in kmin:kmax)
  {
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

      if(subnetScore=="mean"){
        sumStat<-c(mean(gList[sample(length(gList$gene),k,replace=TRUE),]$score,na.rm=TRUE))
      }
      if(subnetScore=="sum") {
        sumStat<-c(sum(gList[sample(length(gList$gene),k,replace=TRUE),]$score,na.rm=TRUE))
      }
      ba<-rbind(ba,sumStat)
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
