### Example

library(roxygen2)
library(devtools)
document()

library(signet)
data(keggPathways)
data(zScores)

bkgd <- backgroundDist(pathwaysList = keggPathways,
                       scores = zScores,
                       kmax = 10)
# data(bkgd)


signetObject <- searchSubnet(pathway = keggPathways,
                             scores = zScores,
                             nullDist = bkgd)


results <- correctOverlap(signetObject,
                          cluster = "max",
                          multipleTesting = TRUE,
                          threshold = 0.05)

writeResults(results)


library(graph)
data("MAPKsig")

set.seed(1)
g1 <- randomEGraph(paste(1:70),0.05)
plot(g1)

randomScore<-data.frame(nodes(g1),(rnorm(length(nodes(g1)))))

bkgd<-matrix(NA,ncol=3,nrow=50)
colnames(bkgd)<-c("k","mu","sigma")

for(k in 1:50)
{
  ba<-NULL
  for(i in 1:10000)
  {
    ba<-c(ba,mean((rnorm(k))))
  }
  bkgd[k,]<-c(k,mean(ba),sd(ba))
}
bkgd<-as.data.frame(bkgd)


randomScore[c(41,42,36,50,6,31,45,14,55,66),2] <-
  sort((rnorm(10000)),decreasing=TRUE)[1:10]+5

signetObject <- searchSubnet(pathway = g1,
                             nullDist = bkgd,
                             scores = randomScore,
                             subnetScore = "mean",
                             iterations = 3000,
                             temperature=0.9995,
                             diagnostic=TRUE,animPlot = 5000)

Rprof(NULL)
summaryRprof("file.out")


plotSubnet(g1,signetObject$table)





hist(replicate(1000,sum(sample(runif(1000),10))))
hist(zScores$z)
