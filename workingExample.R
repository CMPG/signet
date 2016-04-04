### Example

library(roxygen2)
library(devtools)
document()

library(signet)
data(keggPathways)
data(zScores)

rm(list=ls())
#select graphs
library(graphite)
reactome<-pathways("hsapiens", "reactome")
kegg<-pathways("hsapiens", "kegg")
nci<-pathways("hsapiens", "nci")

nodes(pathwayGraph(reactome[[1]]))
#choose type of gene identifier
kegg<- convertIdentifiers(kegg, "entrez")
reactome<-convertIdentifiers(reactome, "entrez")
nci<-convertIdentifiers(nci, "entrez")

# convert to graphNEL format
kegg<-lapply(kegg,pathwayGraph)
attr(kegg,"database")<-"kegg"

reactome<-lapply(reactome,pathwayGraph)
attr(reactome,"database")<-"reactome"

nci<-lapply(nci,pathwayGraph)
attr(nci,"database")<-"nci"

allgraphs<-c(kegg,reactome,nci)

# save(kegg,reactome,nci,file="allgraphs.rda")

reactome_summary<-graphSummary(reactome)
kegg_summary<-graphSummary(kegg)
nci_summary<-graphSummary(nci)

save(kegg_summary,reactome_summary,nci_summary,file="pathways_summary_entrez.rda")

kegg_clean<-filterGraphs(kegg,kegg_summary,10,0.5)
reactome_clean<-filterGraphs(reactome,reactome_summary,10,1)
nci_clean<-filterGraphs(nci,nci_summary,10,0.5)
length(kegg_clean)

save(kegg_clean,reactome_clean,nci_clean,file="allgraphs_clean_entrez.rda")

bkgd <- backgroundDist(pathwaysList = keggPathways,
                       scores = zScores,
                       subnetScore = "mean",
                       kmax = 50,
                       iterations = 200)

# save(bkgd,file="kegg_bkgd.rda")

signetObject <- searchSubnet(pathway = keggPathways[[2]],
                             scores = zScores,
                             subnetScore = "mean",
                             iterations = 2000,
                             replicate = 1,
                             nullDist = bkgd,
                             temperature=0.9995,
                             diagnostic=TRUE,
                             animPlot = 1000)

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
    sc<-mean(rnorm(k))-mean(tail(sort(mean(rnorm(50))),5))
    ba<-c(ba,sc)
  }
  bkgd[k,]<-c(k,mean(ba),sd(ba))
  print(k)
}
bkgd<-as.data.frame(bkgd)


randomScore[c(41,42,36,50,6,31,45,14,55,66,51),2] <- rnorm(11,mean=5)
library(signet)

help(signet)
signetObject <- searchSubnet(pathway = g1,
                             nullDist = bkgd,
                             scores = randomScore,
                             subnetScore = "delta",
                             iterations = 5000,
                             temperature=0.9995,
                             diagnostic=TRUE,
                             animPlot = 5000)

