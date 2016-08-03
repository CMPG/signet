#' Compute a subnetwork score

computeScore <- function(signetObject,score="mean")
{
  if(score == "mean") {
    subnetStat <- mean(signetObject[signetObject$state,]$score)
  }
  if(score == "sum") {
    subnetStat <- sum(signetObject[signetObject$state,]$score)
  }
  if(score == "ideker") {
    k<-length(signetObject[signetObject$state,]$score)
    subnetStat <- (1/sqrt(k))*sum(signetObject[signetObject$state,]$score)
  }
  if(score == "delta") {
    subnetStat <- mean(signetObject[signetObject$state,]$score)-mean(
      tail(sort(signetObject[!signetObject$state,]$score),5))
  }
  return(subnetStat)
}

