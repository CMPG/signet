#' Null distribution
#'
#' Generate the high-scores null distribution to compute empirical p-values
#' for each pathway tested.
#'
#' @param pathwaysList A list of graphNEL objects.
#' @param scores A data frame in which the first column if the gene list
#' and the second column contains the gene scores.
#' @param iterations Number of null high-scores to compute.
#'
#' @keywords simulated annealing, null distribution
#' @export

nullDist<-function(pathwaysList,
                         scores,
                         iterations = 1000,
                         bkgd) {
  requireNamespace("graph",quietly=TRUE)
  colnames(scores) <- c("gene","score")

  genesDB <- unique(unlist(sapply(pathwaysList,graph::nodes)))
  scores <- scores[scores$gene %in% genesDB,]

  if (missing(bkgd)) {
    bkgd <- data.frame(k=1:200,
                             mu=sqrt(1:200)*mean(scores$score,na.rm=TRUE),
                             sigma=rep(sd(scores$score,na.rm=TRUE),200))
  }


  cat("  Generating null distribution...\n")
  bk<-lapply_pb(1:iterations,function(x){

      newscores<-data.frame(gene=scores$gene,score=sample(scores$score))
      cond<-TRUE
      while(cond){
        grap<-sample(1:length(pathwaysList),1)
        HSS<-try(searchSubnet(pathwaysList[[grap]],
                              scores=newscores,
                              bkgd,iterations=1000,verbose=FALSE),silent=TRUE)
        if(class(HSS)!="try-error") cond<-FALSE
      }
      return(HSS$subnet_score)

  })

  cat("\n  Done!\n\n")
  return(unlist(bk))
}
