### Example

# library(roxygen2)
# library(devtools)
# document()

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
help(graph)
data(MAPKsig)
nat = rep(TRUE, length(nodes(MAPKsig)))
names(nat) = nodes(MAPKsig)
plot(MAPKsig)

MAPKsig<-graphNEL(MAPKsig)
g1 = randomEGraph(paste(1:50), edges=120)
plot(g1)

randomScore<-data.frame(nodes(g1),rnorm(length(nodes(g1))))
randomScore[c(1,8,2,6,49,16,15),2] <- sort(rnorm(1000),decreasing = TRUE)[1:7]+1

signetObject <- searchSubnet(pathway = g1,
                             nullDist = data.frame(k=c(1:50),mu=0,sigma=seq(1,0.2,length.out = 50)),
                             scores = randomScore,
                             iterations = 5000,
                             temperature=0.999,
                             directed=FALSE,
                             diagnostic=TRUE)

plotSubnet(MAPKsig,signetObject$table)
