#' Null distribution
#'
#' Generate the high-scores null distribution to compute empirical p-values
#' for each pathway tested.
#'
#' @param pathways A list of graphNEL objects.
#' @param scores A data frame in which the first column is the gene list
#' and the second column contains the gene scores.
#' @param n Number of null high-scores to compute.
#' @param background Background distribution.
#'
#' @keywords simulated annealing, null distribution
#' @export
#' @examples
#' \dontrun{
#' data(daub13)
#'
#' nullDistrib <- nullDist(kegg_human, scores, n = 1000)
#' }

nullDist<-function(pathways,
                   scores,
                   n = 1000,
                   background) {

  requireNamespace("graph",quietly=TRUE)
  colnames(scores) <- c("gene","score")

  if (missing(background)) {

    background <- data.frame(k=1:200,
                             mu=sqrt(1:200)*mean(scores$score,na.rm=TRUE),
                             sigma=rep(sd(scores$score,na.rm=TRUE),200))

  }

  genesDB <- unique(unlist(sapply(pathways,graph::nodes)))
  scores <- scores[scores$gene %in% genesDB,]

  cat("  Generating null distribution...\n")

  bk<-lapply_pb(1:n,function(x){

      newscores<-data.frame(gene=scores$gene,score=sample(scores$score))
      cond<-TRUE

      while(cond){

        grap<-sample(1:length(pathways),1)

        HSS<-try(searchSubnet(pathways[[grap]],
                              scores=newscores,
                              background,
                              iterations=1000,
                              verbose=FALSE),silent=TRUE)

        if(class(HSS)!="try-error") cond<-FALSE

      }
      return(HSS$subnet_score)

  })

  cat("\n  Done!\n\n")

  return(unlist(bk))
}
