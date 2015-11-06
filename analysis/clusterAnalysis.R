

allResRmax<-list()
for(i in 1:1548)
{
  filename<-paste("//cmpgneo.unibe.ch/gouy/datambe/results/Routput",i,".rda",sep="")
  bla<-try(load(filename))
  if(class(bla)!="try-error") allResRmax[[i]]<-output
  else allResRmax[[i]]<-NULL
}
Rnullmax<-unlist(read.table("//cmpgneo.unibe.ch/gouy/datambe/null/Rnull.txt"))
write.table(returnTable(allResRmax,reactome,Rnullmax),file="RMax.tsv",
            sep="\t",quote=FALSE,row.names=FALSE)




allResKmean<-list()
for(i in 1:250)
{
  filename<-paste("//cmpgneo.unibe.ch/gouy/meanFst/results/Koutput",i,".rda",sep="")
  bla<-try(load(filename))
  if(class(bla)!="try-error") allResKmean[[i]]<-output
  else allResKmean[[i]]<-NULL
}
Knull<-unlist(read.table("//cmpgneo.unibe.ch/gouy/meanFst/null/Knull.txt"))
write.table(returnTable(allResKmean,kegg,Knull),file="KeggMean.tsv",
            sep="\t",quote=FALSE,row.names=FALSE)

allResKmax<-list()
for(i in 1:250)
{
  filename<-paste("//cmpgneo.unibe.ch/gouy/maxFst/results/Koutput",i,".rda",sep="")
  bla<-try(load(filename))
  if(class(bla)!="try-error") allResKmax[[i]]<-output
  else allResKmax[[i]]<-NULL
}
Knullmax<-unlist(read.table("//cmpgneo.unibe.ch/gouy/maxFst/null/Knull.txt"))
write.table(returnTable(allResKmax,kegg,Knullmax),file="KeggMax.tsv",
            sep="\t",quote=FALSE,row.names=FALSE)






save(allResR,file="maxFstReactome.rda")
load(file="QUEBEC_FINAL.rda"))


returnTable<-function(outputSignet,pDatabase,nullDistribution)
{
  subnetSize<-unlist(lapply(outputSignet,function(x) {
    stat<-x$size;if(is.null(stat)) stat<-NA; return(stat)
    }))
  netSize<-unlist(lapply(outputSignet,function(x) {
    stat<-length(x$table$gene);if(is.null(stat)) stat<-NA; return(stat)
    }))
  subnetScore<-unlist(lapply(outputSignet,function(x) {
    stat<-x$score;if(is.null(stat)) stat<-NA; return(stat)
    }))

  pvalues<-unlist(lapply(subnetScore,function(x) {
    stat<-mean(nullDistribution>x,na.rm=TRUE);if(is.null(stat)) stat<-NA; return(stat)
    }))
  pvalues[subnetSize<3]<-NA

  pathwaysNames<-names(pDatabase)

  require(qvalue)
  qvalues<-qvalue(pvalues)$qvalues

  geneListEntrez<-unlist(lapply(outputSignet,function(x) {
    stat<-paste(x$table[x$table$state,]$gene,collapse=";");
    if(is.null(stat)) stat<-NA; return(stat)
  }))
  url<-unlist(lapply(outputSignet,function(x) {
    stat<-paste("http://www.ncbi.nlm.nih.gov/gene/?term=",
                paste(x$table[x$table$state,]$gene,collapse="+"),collapse="");
    if(is.null(stat)) stat<-NA; return(stat)
  }))

  out<-data.frame(pathwaysNames,netSize,subnetSize,subnetScore,
                  pvalues,qvalues,geneListEntrez,url)


  return(out)
}

write.table(returnTable(allResRmax,reactome,Rnullmax),file="new.tsv",
            sep="\t",quote=FALSE,row.names=FALSE)


subnetScore<-unlist(lapply(allResRmax,function(x) {
  stat<-x$score;if(is.null(stat)) stat<-NA; return(stat)
}))
subnetSize<-unlist(lapply(allResRmax,function(x) {
  stat<-x$size;if(is.null(stat)) stat<-NA; return(stat)
}))
geneListUnique<-unlist(lapply(allResRmax[subnetScore>7&subnetSize>2],function(x) {
  stat<-paste(x$table[x$table$state,]$gene);
  if(is.null(stat)) stat<-NA; return(stat)
}))


Rnullmax<-unlist(read.table("//cmpgneo.unibe.ch/gouy/maxFst/null/AllRnull.txt"))
pvalues<-unlist(lapply(subnetScore,function(x) {
  stat<-mean(Rnullmax>x,na.rm=TRUE);if(is.null(stat)) stat<-NA; return(stat)
}))
pvalues[subnetSize<3]<-NA
hist(pvalues)
pathwaysNames<-names(pDatabase)

require(qvalue)
qvalues<-qvalue(pvalues)
max(qvalues$pvalues[qvalues$qvalues<0.05],na.rm=TRUE)
hist(qvalue(pvalues))

#cluster
overlapReact<-matrix(0,nrow=length(allResRmax),ncol=length(allResRmax))
for(i in 1:length(allResRmax))
{
  for(j in 1:length(allResRmax))
  {
    gi<-allResRmax[[i]]$table[allResRmax[[i]]$table$state,]$gene
    gj<-allResRmax[[j]]$table[allResRmax[[j]]$table$state,]$gene
    inter<-length(intersect(gi,gj))
    un<-length(union(gi,gj))
    if(!is.null(inter)&!is.null(un)) overlapReact[i,j]<-inter/un
    else overlapReact[i,j]<-NA
  }
  print(i)
}

# head(overlapReact)
# torem<-which(colSums(overlapReact)<c(-500) & rowSums(overlapReact)<c(-500))
# torem2<-which(!is.na(subnetSize) & subnetSize > 2 & subnetScore > 8.3)

torem2<- !is.na(pvalues) & pvalues < 0.05
# overlap2<-overlapReact[-c(torem,torem2),-c(torem,torem2)]
overlap2<-overlapReact[c(torem2),c(torem2)]


# overlapReact[is.nan(overlapReact)]<- -1
# overlapReact[is.na(overlapReact)]<- -1

dim(overlap2)
d<-as.dist(100*(1-overlap2))

library(graphite)
names(reactome)[torem2]

h = hclust(d)
plot(h, main = "Overlapping percentage between subnetworks", sub = "", xlab="",
     axes = FALSE, hang = -1,labels=names(reactome)[torem2],cex=0.7)
lines(x = c(0,0), y = c(0,100), type = "n") # force extension of y axis
axis(side = 2, at = seq(0,100,10), labels = seq(100,0,-10))
abline(h=c(70,50),col=2,lty=2,lwd=2)

plot(hclust(d),labels=FALSE)
?hclust

nclust<-NULL
for(i in 1:100)
{
  groups<-cutree(h,h=100-i)
  nclust<-c(nclust,max(groups))
}

plot(1:100,nclust)

pvalues[!is.na(pvalues)][o]
o<-order(pvalues[!is.na(pvalues)])
cutree(h,h=80)[o]

qv<-qvalue(tapply(pvalues[!is.na(pvalues)],as.factor(cutree(h,h=50)),max))
hist(qv)
max(qv$pvalues[qv$qvalues<0.05])
unlist(qv$qvalues)

newRes<-allResRmax[torem2]

sco<-unlist(lapply(newRes,function(x) return(x$score)))

fina<-unlist(tapply(sco,as.factor(groups),function(x) c(which(sco==max(x)))))

unlist(lapply(newRes[fina],function(x) return(x$score)))





GGlist<-lapply(test,strsplit,split=",")

overlapExp<-matrix(0,nrow=41,ncol=41)
for(i in 1:41)
{
  for(j in 1:41)
  {
    gi<-GGlist[[i]][[1]]
    gj<-GGlist[[j]][[1]]
    inter<-length(intersect(gi,gj))
    un<-length(union(gi,gj))
    if(!is.null(inter)&!is.null(un)) overlapExp[i,j]<-inter/un
    else overlapExp[i,j]<-NA
  }
  print(i)
}


overlapExp<-overlapReact
d<-as.dist(100*(1-overlapExp))

library(graphite)
tnames<-readClipboard()

h = hclust(d)
plot(h, main = "Overlapping percentage between tissues", sub = "", xlab="",
     axes = FALSE, labels=tnames,hang = -1,cex=0.7)
lines(x = c(0,0), y = c(0,100), type = "n") # force extension of y axis
axis(side = 2, at = seq(0,100,10), labels = seq(100,0,-10))








diag(overlapReact)<-NA
getOverlap<-function(Graph)
{
  overlap<-matrix(0,nrow=length(Graph),ncol=length(Graph))
  for(i in 1:length(Graph)){for(j in 1:length(Graph)){
    ov1<-sum(Graph[[j]]@nodes %in% Graph[[i]]@nodes)/length(c(Graph[[i]]@nodes))
    ov2<-sum(Graph[[i]]@nodes %in% Graph[[j]]@nodes)/length(c(Graph[[j]]@nodes))
    overlap[i,j]<-max(ov1,ov2)
  }}
  prop<-NULL
  for(i in seq(0,1,0.01))
  {
    prop<-c(prop,mean(overlap[lower.tri(overlap)]>=i,na.rm=TRUE))
  }
  return(prop)
}



source("https://bioconductor.org/biocLite.R")
library("reactome.db")

gl<-readClipboard()

xx <- as.list(reactomeEXTID2PATHID)

reactid<-unique(unlist(xx[gl[!is.na(gl)]]))
#to GO
x2 <- as.list(reactomeREACTOMEID2GO)
go<-unique(unlist(x2[reactid]))



geneout<-c(unique(geneListUnique))

bool<-as.character(annot[annot$entrezID %in% geneout,]$query)
subse<-(dat[dat$ensembleID %in% bool,])
head(subse)
sco<-unlist(tapply(subse$score,subse$ensembleID,max))

plot(subse[subse$score%in%sco,]$freq.core,subse[subse$score%in%sco,]$freq.front)
hist(subse$score)

subse[subse$score%in%sco,]$
head(subse)
colnames(annot)


paste("http://www.ensembl.org/Homo_sapiens/Location/View?db=core;r=",subse[subse$score%in%sco,]$chrom.num,":",subse[subse$score%in%sco,]$pos,";vdb=variation",sep="")

length(allRes)

subnetSize<-unlist(lapply(allResRmean,function(x) {stat<-x$size;if(is.null(stat)) stat<-NA; return(stat)}))
netSize<-unlist(lapply(allResR,function(x) {stat<-length(x$table$gene);if(is.null(stat)) stat<-NA; return(stat)}))
subnetScore<-unlist(lapply(allResRmean,function(x) {stat<-x$score;if(is.null(stat)) stat<-NA; return(stat)}))

hist(netSize)
hist(subnetSize)
hist(subnetScore,freq = FALSE)
plot(subnetScore,subnetSize,pch=16,cex=0.5)

url<-unlist(lapply(allResRmean,function(x) {
  stat<-paste("http://www.ncbi.nlm.nih.gov/gene/?term=",
              paste(x$table[x$table$state,]$gene,collapse="+"),collapse="");
  if(is.null(stat)) stat<-NA; return(stat)
}))

(url[which(subnetSize>2 & pvalues<0.05)])

###test
Rnull<-unlist(read.table("//cmpgneo.unibe.ch/gouy/meanFst/null/Rnull.txt"))
plot(density(Rnull))
lines(density(Knull))
abline(v=subnetScore,col=rgb(0,0,0,alpha=0.1))
abline(v=quantile(Rnull,probs=0.95),col=2,lwd=2)
abline(v=quantile(Knull,probs=0.999),col=2,lwd=2)

pvalues<-unlist(lapply(subnetScore,function(x) {stat<-mean(Rnull>x,na.rm=TRUE);if(is.null(stat)) stat<-NA; return(stat)}))

hist(pvalues[!is.na(pvalues)],breaks=10)

length(subnetScore[!is.na(subnetScore)])

names(kegg[which(pvalues<0.05)])
#before multiple testing correction: correct for overlapping

##overlapping proportion in the results
#percentage of genes appearing more than once
mean(table(unlist(lapply(allRes[which(pvalues<0.01)],function(x)
  return(x$table[x$table$state,]$gene))))>1)

cat(names(table(unlist(lapply(allRes[which(pvalues<0.01)],function(x)
  return(x$table[x$table$state,]$gene))))))

###correction

require(qvalue)
qvalues<-qvalue(pvalues[subnetSize>2&!is.na(pvalues)])
hist(qvalues)
candidates<-which(qvalues$qvalues<0.05 & pvalues[subnetSize>2&!is.na(pvalues)]<0.01)
names(allGraph[pvalues<0.01&subnetSize>2])[!is.na(names(allGraph[pvalues<0.01&subnetSize>2]))]
which(pvalues<0.01&subnetSize>2)[16]

names(reactome[which(pvalues<0.05)])

allRes[[857]]




hist(subnetScore[1549:1800])
allRes[[1548+157]]

nullDistribut<-unlist(read.table("nullDistribution.txt"))
hist(nullDistribut)

length(nullDistribut)

#check for uniformity
pvalnull<-array(NA,400)
for(i in 1:400)
{
  pvalnull[i]<-mean(nullDistribut[!is.na(nullDistribut)][i]<nullDistribut[!is.na(nullDistribut)][])
}
hist(pvalnull)



hist(subnetScore[subnetSize>2&!is.na(pvalues)],freq = FALSE,ylim=c(0,0.2))
plot(density(nullDistribut)$x,density(nullDistribut)$y,col=1,lwd=2,
     lty=2,type="l",ylab="Density",xlab="Subnetwork score")

lines(density(subnetScore[subnetSize>2&!is.na(pvalues)],na.rm=TRUE)$x,density(subnetScore[subnetSize>2&!is.na(pvalues)],na.rm=TRUE)$y,col=3,lwd=2)
legend("topright",lty=c(2,1),col=c(1,3),lwd=2,c("Null","Observed"))

