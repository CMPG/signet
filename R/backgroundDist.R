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
#' @param iterations Number of iterations to make the background distribution.
#'
#' @keywords simulated annealing, background distribution
#' @export

backgroundDist<-function(pathwaysList,
                         scores,
                         kmin = 1,
                         kmax = 100,
                         iterations = 1000) {
  requireNamespace("graph",quietly=TRUE)
  if(missing(pathwaysList)) {
    pathwaysList <- NULL
    gList <- scores
  }

  colnames(scores) <- c("gene","score")

  if(kmax=="max") { # kmax not specified: max k in the graph list
    kmax <- max(unlist(lapply(pathwaysList,function(x)length(nodes(x)))))
  }

  cat("  Generating background distribution...\n")

  bk<-lapply_pb(kmin:kmax,function(x){
    ba<-NULL

    ba<-unlist(lapply(1:iterations,function(y){

      if(length(pathwaysList)>0) {

        glis<-NULL
        while(length(glis)<x+5)
        {
          path<-pathwaysList[[sample(length(pathwaysList),1)]]
          glis<-graph::nodes(path)
        }
        gList<-scores[scores$gene %in% glis,]
      }

      gList2<-gList[sample(length(gList$gene),x,replace=TRUE),]
      sumStat <- (1/sqrt(x))*sum(gList2$score)

      return(sumStat)
    }))

    return(c(k=x,mu=mean(ba,na.rm=TRUE),sigma=sd(ba,na.rm=TRUE)))
  })

  nullD<-as.data.frame(do.call("rbind",bk))
  cat("\n  Done!\n\n")
  return(nullD)
}
