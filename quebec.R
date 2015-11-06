source("https://bioconductor.org/biocLite.R")
biocLite("mygene")
library("mygene")#problems with dependencies versions... -> reinstall all BioC packages

dat<-read.table("./analysis/quebec.txt",header=TRUE,sep="\t")
dat$ensembleID<-substring(as.character(dat$ens.gene),1,15)
#ensemble ID -> entrez ID +
#
# out <- queryMany(dat$ensembleID, scopes="ensembl.gene", fields="symbol", species="human")
# out2 <- queryMany(dat$ensembleID, scopes="ensembl.gene", fields=c("entrezgene","name","alias","summary"), species="human")
#
# alias<-lapply(out2$alias,function(x) paste(x, collapse = ";"))
# alias<-unlist(alias)
#
# out2$alias<-alias
# out3<-out2[!duplicated(out2$entrezgene), ]
# dim(out3)
#
# write.table(out3,file="gene_annotation.tsv", sep="\t", row.names = FALSE, quote = FALSE)
#
#
# duplist<-names(table(out3$query)[table(out3$query)>1])
# #130 genes have more than one entrez ID
# out4<-out3[out3$query%in%duplist,]
# out5<-NULL
# for(gene in 1:length(duplist))
# {
# #   print(duplist[gene])
#   out5<-rbind(out5,paste(as.data.frame(out4[out4$query==duplist[gene],])))
# }
#
# colnames(out5)<-colnames(out6)
#
# out6<-out3[!out3$query%in%duplist,]
# out6<-rbind(out6,out5)
#
# write.table(out6,file="gene_annotation_unique.tsv", sep="\t", row.names = FALSE, quote = FALSE)



annot<-read.table("./analysis/gene_annotation_unique.tsv", sep="\t",header=TRUE,na.strings="NA",quote = '"')
annot<-annot[annot$keep==1,]
annot$entrezID<-as.numeric(as.character(annot$entrezID))

entrezID<-lapply(dat$ensembleID,function(x) annot[annot$query==x,][1,]$entrezgene)
# name<-lapply(dat$ensembleID,function(x) annot[annot$query==x,][1,]$name)
# summa<-lapply(dat$ensembleID,function(x) annot[annot$query==x,][1,]$summary)

dat$entrezID <- entrezID
# dat$symbol <- out2$symbol
# dat$name <- name
# # dat$alias <- alias
# dat$summary <- summa

## compute delta freq
within<-((1-(dat$freq.front/100)^2)+(1-(dat$freq.core/100)^2))/2
between<-(1-((dat$freq.front-dat$freq.front)/200)^2)
fst<-1-(within/between)
hist(fst)
dat$score<-fst#c(dat$freq.front-dat$freq.core)

# save(dat,file="quebecRawData.rda")
load(file="quebecRawData.rda")

###tronche des stats
plot(dat$score~as.factor(dat$annovar))


plot(dat$gerp,dat$score,cex=0.5,pch=16,col=rgb(.0,.0,.0,alpha=0.1))

snpPerGene<-tapply(dat$score,dat$ensembleID,length)
maxScore<-tapply(dat$score,dat$ensembleID,max)
meanScore<-tapply(dat$score,dat$ensembleID,mean,na.rm=TRUE)
q90Score<-tapply(dat$score,dat$ensembleID,function(x) mean(x[x>quantile(x,probs=c(0.9),na.rm=TRUE)]))
quant80Score<-tapply(dat$score,dat$ensembleID,function(x) quantile(x,probs=c(0.8),na.rm=TRUE))
quant90Score<-tapply(dat$score,dat$ensembleID,function(x) quantile(x,probs=c(0.9),na.rm=TRUE))

par(mfrow=c(2,2))
hist(snpPerGene,breaks=50)
hist(maxScore)
hist(meanScore)
hist(quant80Score)


plot(log10(snpPerGene),maxScore,cex=0.5,pch=16,col=rgb(.0,.0,.0,alpha=0.1),
     main=format(cor(log10(snpPerGene),maxScore),digits=3))
plot(log10(snpPerGene),meanScore,cex=0.5,pch=16,col=rgb(.0,.0,.0,alpha=0.1),
     main=format(cor(log10(snpPerGene),meanScore),digits=3))
plot(log10(snpPerGene),q90Score,cex=0.5,pch=16,col=rgb(.0,.0,.0,alpha=0.1),
     main=format(cor.test(log10(snpPerGene),q90Score)$estimate,digits=3))
plot(log10(snpPerGene),quant80Score,cex=0.5,pch=16,col=rgb(.0,.0,.0,alpha=0.1),
     main=format(cor.test(log10(snpPerGene),quant80Score)$estimate,digits=3))


plot(meanScore,maxScore,cex=0.5,pch=16,col=rgb(.0,.0,.0,alpha=0.1))
plot(meanScore,q90Score,cex=0.5,pch=16,col=rgb(.0,.0,.0,alpha=0.1))
plot(maxScore,q90Score,cex=0.5,pch=16,col=rgb(.0,.0,.0,alpha=0.1))


##get scores for the gene list
bla<-NULL
gene<-lapply(names(quant90Score),function(x){bla<-c(annot[annot$query==x,]$entrezID);
    if(sum(bla,na.rm = TRUE)<1) bla<- NA;
  # if(is.na(bla)) bla<- c(-1);
  return(bla)})

gene<-unlist(gene)
scores<-data.frame(gene=gene,score=quant90Score)

save(scores,file="./analysis/fstQuebec.rda")
hist(scores$score)

load("./analysis/scoresQuebec.rda")


#now we can run the analysis
library(signet);library(devtools)
document()
data("keggPathways")

quebecNull2<-nullDistribution(keggGraph,scores = scores,kmin=2,kmax=4,iterations=10000)
for(i in 1:450) quebecNullKegg<-rbind(quebecNullKegg,quebecNullKegg[49,])
quebecNullKegg$k<-2:500

save(quebecNullKegg,file="./analysis/quebecBackgroundKegg.rda")

load("quebecNull.rda")
#before generating the entire null distribution

# save(resultsQuebec,file="resultsQuebec.rda") #WITH DELTA F
# load(file="resultsQuebec.rda") #WITH DELTA F

all3<-searchSubnet(keggGraph[101:length(keggGraph)],
                 scores=scores,
                 nullDist = quebecNullKegg,
                 iterations = 4000,
                 replicates = 2,
                 temperature = 0.999,
                 diagnostic=TRUE,
                 verbose=FALSE)

allKEGG<-c(all[1:9],list(NULL),all[10:99],all3)

#####
### Test the significance of the subnetworks
#####

# compute the null distribution. Permute the scores and apply the algorithm N times
iter<-1000
# test<-array(NA,iter)
#plot null distribution
hist(test,
     main="Null distribution",
     xlab="Subnetwork score",
     breaks=5)


#check for uniformity
pvalnull<-array(NA,length(test[!is.na(test)]))
for(i in 1:150)
{
  pvalnull[i]<-mean(test[!is.na(test)][i]<test[!is.na(test)])
}
hist(pvalnull)

# compute empirical p-value and return the pathways IDs for which FDR < 0.001
subnetSize<-unlist(lapply(allKEGG,function(x) {stat<-x$size;if(is.null(stat)) stat<-NA; return(stat)}))
netSize<-unlist(lapply(allKEGG,function(x) {stat<-length(x$table$gene);if(is.null(stat)) stat<-NA; return(stat)}))
subnetScore<-unlist(lapply(allKEGG,function(x) {stat<-x$score;if(is.null(stat)) stat<-NA; return(stat)}))
hist(subnetScore)

pvalues<-unlist(lapply(subnetScore,function(x) {stat<-mean(test>x,na.rm=TRUE);if(is.null(stat)) stat<-NA; return(stat)}))
hist(pvalues)

require(qvalue)
qvalues<-qvalue(pvalues,fdr.level=0.05)$qvalues

candidates<-which(qvalue(pvalues)$qvalues<0.1)
names(kegg[candidates])

resultsKEGG<-data.frame(name=names(kegg),qvalues=qvalue(pvalues)$qvalues,pvalue=pvalues,subnetScore=subnetScore,
           subnetSize=subnetSize,netSize=netSize)

write.table(resultsKEGG,quote = FALSE,row.names = FALSE,file="quebecKEGG.tsv",sep="\t")


save(allKEGG,resultsKEGG,test,quebecNullKegg,file="auebecKEGGall.rda")

##overlapping proportion in the results
#percentage of genes appearing more than once
mean(table(unlist(lapply(allKEGG[which(qvalue(pvalues)$qvalues<0.05)],function(x)
  return(x$table[x$table$state,]$gene))))>3)





#####
### Procedure for pathways overlapping correction
#####

library(graphite)

react<-pathways("hsapiens","reactome")
react<-convertIdentifiers(react, "entrez")
reactGraph<-lapply(react,pathwayGraph,edge.types=NULL)

kegg<-pathways("hsapiens","kegg")
kegg<-convertIdentifiers(kegg, "entrez")
keggGraph<-lapply(kegg,pathwayGraph,edge.types=NULL)

biocarta<-pathways("hsapiens","biocarta")
biocarta<-convertIdentifiers(biocarta, "entrez")
biocartaGraph<-lapply(biocarta,pathwayGraph,edge.types=NULL)

panther<-pathways("hsapiens","panther")
panther<-convertIdentifiers(panther, "entrez")
pantherGraph<-lapply(panther,pathwayGraph,edge.types=NULL)

nci<-pathways("hsapiens","nci")
nci<-convertIdentifiers(nci, "entrez")
nciGraph<-lapply(nci,pathwayGraph,edge.types=NULL)


all<-searchSubnet(nciGraph,
                  scores=scores,
                  nullDist = quebecNull,
                  iterations = 5000,
                  replicates = 3,
                  temperature = 0.9995,
                  diagnostic=TRUE,
                  verbose=FALSE)


#pathway object to graphNEL

overlap<-matrix(0,nrow=length(reactGraph),ncol=length(reactGraph))
#we build a matrix with the pairwise overlap between pathways
for(i in 1:length(reactGraph))
{
  for(j in 1:length(reactGraph))
  {
    overlap[i,j]<-sum(reactGraph[[i]]@nodes %in% reactGraph[[j]]@nodes)/length(c(reactGraph[[j]]@nodes))
  }
  print(i)
}

hist(overlap[overlap>0.0])

diag(overlap)<-0

sum(overlap>0.5,na.rm=TRUE)

mat<-matrix(NA,ncol=6,nrow=length(reactGraph)^2)
for(i in 1:length(reactGraph))
{
  for(j in 1:length(reactGraph))
  {
    line<-c(i,j,length(reactGraph[[i]]@nodes),length(reactGraph[[j]]@nodes),
            c(overlap[i,j]),
            sum(reactGraph[[i]]@nodes %in% reactGraph[[j]]@nodes))
    mat[i+(j-1)+(i-1)*length(reactGraph),]<-line
  }
  print(i)
}

mat2<-mat[mat[,5]>0.1&!is.na(mat[,5]),]

colnames(mat2)<-c("pathI","pathJ","sizeI","sizeJ","overlap","commonGenes")

save(mat,mat2,react,reactGraph,overlap,file="overlapReactome.rda")
load(file="overlapReactome.rda")

mat3<-as.data.frame(mat2)
#case1
case<-mat3[mat3$overlap==1,]

case1<-case[case$pathJ%in%as.numeric(names(table(case$pathJ)[table(case$pathJ)>1])),]

toremove<-as.numeric(unique(case1$pathI))
#case 1 and 2: juste have to remove the p. we dont care about the names

case3<-case[case$pathJ%in%as.numeric(names(table(case$pathJ)[table(case$pathJ)>1])),]


react[[131]]
superpathways
toremove<-












library(graph)

hierarchy<-read.table("reactomeHierarchy.txt",sep="\t")
hierarchy[,1]<-as.character(hierarchy[,1])
hierarchy[,2]<-as.character(hierarchy[,2])

edL<-list()
for(path in unique(c(hierarchy[,1]))){
  edL[[path]]<-hierarchy[hierarchy[,1]==path,2]
}
x<-graphNEL(nodes=unique(unlist(hierarchy)),edgeL = edL,edgemode = "directed")

leaves<-leaves(x,degree.dir = "out")
parents<-unlist(lapply(leaves,function(x) return(hierarchy[hierarchy[,2]==x,1][1])))

filter<-data.frame(leaves=leaves,parents=parents)

idreactome<-read.table("idreactome.txt",sep="\t")

name<-unlist(lapply(as.character(filter$leaves),
                    function(x) return(as.character(idreactome$V2)[as.character(idreactome$V1)==x][1])))
filter<-cbind(filter,name)
head(filter)

library(graphite)
reactomedb<-pathways("hsapiens","reactome")

for(i in 1:10)
{
  p<-reactomedb[[as.character(filter$name[i])]]
  if(!is.null(p)) x<-(length(nodes(p)))
  else x<-NA
}


ids<-unlist(lapply(reactomedb,function(x) return(x@id)))
length<-unlist(lapply(reactomedb,function(x) return(length(nodes(x)))))

sizes<-unlist(lapply(as.character(filter$leaves),
                    function(x) return(as.numeric(unlist(length))[as.character(unlist(ids))==x][1])))

filter<-cbind(filter,sizes)

tapply(as.character(filter$name),filter$parents,function(x){
  plist<-as.character(x)
  print(plist)
  glist<-NULL
  for(path in plist)
  {
    p<-reactomedb[["path"]]
    if(!is.null(p)) glist<-c(glist,nodes(reactomedb[["path"]]))
    print(p)
    print(mean(table(glist)>1))
  }
})

str(reactomedb[1]@id)
