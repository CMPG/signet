#' Compute a subnetwork score
#' @export

computeScore <- function(signetObject, score = "ideker")
{
  activeNet <- signetObject$network[signetObject$network$active, ]

  if (score == "mean") {
    subnetStat <- mean(activeNet$score)
  }

  if (score == "sum") {
    subnetStat <- sum(activeNet$score)
  }

  if (score == "ideker") {
    k <- length(activeNet$score)
    subnetStat <- (1/sqrt(k))*sum(activeNet$score)
  }

  if (score == "delta") {
    subnetStat <- mean(activeNet$score)-mean(tail(sort(activeNet$score),5))
  }

  return(subnetStat)
}

