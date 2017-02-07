## simulated annealing in c++
library(Rcpp)
cppFunction('NumericVector rowSumsC(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector out(nrow);

  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      total += x(i, j);
    }
    out[i] = total;
  }
  return out;
}')
set.seed(1014)
x <- matrix(sample(100), 10)
rowSums(x)
#>  [1] 458 558 488 458 536 537 488 491 508 528
rowSumsC(x)

Lis<-createSignetObject(kegg_clean[[1]],testWmax[,c(1,6)],100,2)
Adj<-getAdjacencyMatrix(kegg_clean[[1]])

Lis$network$active[1]<-

bkgd<-data.frame(k=1:100,mu=0.2,sigma=1)

Adj

cppFunction('float simAnnCpp(NumericMatrix matNet, NumericMatrix matSA, NumericMatrix matAdj, NumericMatrix matBK, int iter,
StringVector genes) {
  float t = matNet(1,1);

  Rcout << genes(0);
  Rcout << matNet;

  for(int i = 0; i < 1; ++i) {
    Rcout << matAdj << std::endl;
  }

  return t;
}')
simAnnCpp(as.matrix(Lis$network),as.matrix(Lis$simulated_annealing),Adj,as.matrix(bkgd),100,colnames(Adj))



list(Lis)<-list

##inputs: 4 matrices
# signetObject$network, table of scores, genes, active
# adjacency matrix
# bkgd distrib
# simul annealing table, size and score evolution

# n iterations

for (i in 1:iterations) {

  activeNet <- as.character(signetObject$network[signetObject$network$active,]$gene)
  adjSubgraph <- adjMatrix[activeNet,]

  if (length(dim(adjSubgraph))>0) {
    boundaries <- adjSubgraph[apply(adjSubgraph,1,sum)>0,]
    boundaries <- colnames(boundaries[,apply(boundaries,2,sum)>0])
    #remove boundaries already active
    boundaries <- boundaries[!boundaries %in% signetObject$network[signetObject$network$active,]$gene]
  } else {
    boundaries <- NULL
  }

  bla <- 0
  while (bla != 1) {
    ge <- sample(activeNet,1)
    test <- graph::subGraph(as.character(activeNet[activeNet != ge]),signetObject$connected_comp)

    bla <- length(graph::connComp(test))
    if (bla == 1) neighbours<-c(ge)
  }

  probneighbour <- 0
  sampNeighbour <- FALSE
  if (length(activeNet) > kmin) {
    probneighbour <- length(activeNet)/(length(boundaries)+length(activeNet))
    if (probneighbour > runif(1)) {
      sampNeighbour <- TRUE
    }
  }

  if (!maximean & verbose) {
    if (length(activeNet) >= max(nullDist$k)) stop(paste("No null distribution
                                                         for subnetwork
                                                         of size > ",max(nullDist$k)))
  }

  if (length(activeNet) > 1 & length(unique(c(neighbours,boundaries))) > 0) {

    # Sample a new gene to toggle its state
    if (!sampNeighbour) newG <- as.character(sample(unique(c(boundaries)),1))
    else newG <- neighbours
    signetObject$network[which(signetObject$network$gene==newG),]$active <-
      !signetObject$network[which(signetObject$network$gene==newG),]$active

    # Compute the subnet score
    sumStat <- computeScore(signetObject, score = subnetScore)

    # Scale the subnet score
    if (maximean) {
      s2 <- sumStat
    } else {
      s2 <- (sumStat-nullDist[nullDist$k == length(activeNet),]$mu)/nullDist[nullDist$k==length(activeNet),]$sigma
    }

    # Keep or not the toggled gene, acceptance probability
    if (s2 < s) {
      prob <- exp((s2-s)/signetObject$simulated_annealing$temperature[i])
      if (prob > runif(1)) {
        s <- s2
        signetObject$simulated_annealing$score_evolution[i] <- s
        signetObject$simulated_annealing$size_evolution[i] <- sum(
          signetObject$network$active
        )
      } else {
        signetObject$network[which(signetObject$network$gene==newG),]$active <-
          !signetObject$network[which(signetObject$network$gene==newG),]$active
        signetObject$simulated_annealing$score_evolution[i] <- s
        signetObject$simulated_annealing$size_evolution[i] <- sum(
          signetObject$network$active
        )
      }
    } else if (s2 > s) {
      s <- s2
      signetObject$simulated_annealing$score_evolution[i] <- s
      signetObject$simulated_annealing$size_evolution[i] <- sum(
        signetObject$network$active
      )
    }
  }
}




