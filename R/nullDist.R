#' Null distribution
#'
#' Generate the high-scores null distribution to compute empirical p-values
#' for each pathway tested.
#'
#' @param pathwaysList A list of graphNEL objects.
#' @param scores A data frame in which the first column if the gene list
#' and the second column contains the gene scores.
#' @param iterations Number of iterations to make the null distribution.
#'
#' @keywords simulated annealing, null distribution
#' @export

nullDist<-function(pathwaysList,
                         scores,
                         iterations = 1000,
                         bkgd) {
  requireNamespace("graph",quietly=TRUE)
  colnames(scores) <- c("gene","score")

  cat("  Generating null distribution...\n")
  bk<-lapply_pb(1:iterations,function(x){

      newscores<-data.frame(gene=scores$gene,score=sample(scores$score))
      grap<-sample(1:length(pathwaysList),1)

      HSS<-try(searchSubnet(pathwaysList[[grap]],
                            scores=newscores,
                            bkgd))
      if(class(HSS)!="try-error") return(HSS$subnet_score)
      else return(NA)

  })

  cat("\n  Done!\n\n")
  return(bk)
}
