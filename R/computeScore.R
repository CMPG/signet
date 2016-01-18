#' Compute a subnetwork score

computeScore <- function(signetObject,score="mean")
{
  if(subnetScore == "mean") {
    subnetStat <- mean(signetObject[signetObject$state,]$score)
  }
  if(subnetScore == "sum") {
    subnetStat <- sum(signetObject[signetObject$state,]$score)
  }
  if(subnetScore == "delta") {
    subnetStat <- mean(signetObject[signetObject$state,]$score)-mean(
      tail(sort(signetObject[!signetObject$state,]$score),5))
  }
  return(subnetStat)
}

