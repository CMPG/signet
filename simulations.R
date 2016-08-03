### Example

library(roxygen2)
library(devtools)
document()

library(signet)

data(keggPathways)
data(zScores)

library(graph);library(igraph)
data("MAPKsig")

set.seed(1)
g1 <- randomEGraph(paste(1:70),0.05)
plot(g1)

randomScore<-data.frame(genes=nodes(g1),scores=(rnorm(length(nodes(g1)))))

bkgd<-matrix(NA,ncol=3,nrow=70)
colnames(bkgd)<-c("k","mu","sigma")

for(k in 1:70)
{
  ba<-NULL
  for(i in 1:5000)
  {
    sc<-(1/sqrt(k))*sum(rnorm(k))#-mean(tail(sort(mean(rnorm(70))),5))
    ba<-c(ba,sc)
  }
  bkgd[k,]<-c(k,mean(ba),sd(ba))
  print(k)
}
bkgd<-as.data.frame(bkgd)


randomScore[c(41,42,36,50,6,31,45,14,55,66,51),2] <- rnorm(11,mean=2)
library(signet)

help(signet)
signetObject <- searchSubnet(pathway = g1,
                             nullDist = bkgd,
                             scores = randomScore,
                             subnetScore = "ideker",
                             iterations = 100000,
                             temperature=0.99995,
                             diagnostic=TRUE,
                             animPlot = 0)

createSignetObject(randomScore,1000)

# faire varier iterations

# faire varier taille réseau

# faire varier taille sous réseau

# faire varier scores sous réseau en écart types
