source("https://bioconductor.org/biocLite.R")
biocLite("mygene")
library("mygene")#problems with dependencies versions... -> reinstall all BioC packages

dat<-read.table("quebec.txt",header=TRUE,sep="\t")
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



annot<-read.table("gene_annotation_unique.tsv", sep="\t",header=TRUE,na.strings="NA",quote = '"')
annot<-annot[annot$keep==1,]
annot$entrezID<-as.numeric(as.character(annot$entrezID))

head(annot)

entrezID<-lapply(dat$ensembleID,function(x) annot[annot$query==x,][1,]$entrezgene)
name<-lapply(dat$ensembleID,function(x) annot[annot$query==x,][1,]$name)
name<-lapply(dat$ensembleID,function(x) annot[annot$query==x,][1,]$summary)

dat$entrezID <- entrezID
# dat$symbol <- out2$symbol
dat$name <- name
dat$alias <- alias
dat$summary <- summa

## compute delta freq
dat$score<-c(dat$freq.front-dat$freq.core)

maxScore<-tapply(dat$score,dat$ensembleID,max)
bla<-NULL
gene<-lapply(names(maxScore),function(x){bla<-c(annot[annot$query==x,]$entrezID);
  print(bla);
  if(sum(bla,na.rm = TRUE)<1) bla<- NA;
  # if(is.na(bla)) bla<- c(-1);
  return(bla)})


gene<-unlist(gene)

scores<-data.frame(gene=gene,score=maxScore)

save(scores,file="scoresQuebec.rda")
head(scores)



#now we can run the analysis
library(signet);library(devtools)
document()
data("keggPathways")

quebecNullR<-nullDistribution(reactGraph,scores = scores,kmin=2,kmax=50,iterations=100)
save(quebecNull,file="quebecNull.rda")
load("quebecNull.rda")
#before generating the entire null distribution
for(i in 1:450) quebecNullR<-rbind(quebecNullR,quebecNullR[49,])
quebecNullR$k<-2:500
plot(quebecNullR)

quebecKegg<-lapply(reactGraph[1:10],searchSubnet,
                         scores=scores,
                         nullDist = quebecNullR,
                         iterations = 5000,
                         replicates = 5,
                         temperature = 0.999,
                         diagnostic=TRUE)
#bug with kegg pathway 106

resultsQuebec<-list()
resultsQuebec<-c(quebecKegg,quebecKegg2,quebecKegg3b,
            quebecKegg3a1,quebecKegg3a2,quebecKegg3a3,
            quebecKegg3c,quebecKegg3d,quebecKegg3e,quebecKegg4,quebecKegg5)
save(resultsQuebec,file="resultsQuebec.rda")
load(file="resultsQuebec.rda")
quebecKegg3a1<-lapply(keggPathways[101:105],searchSubnet,
                   scores=scores,
                   nullDist = quebecNull,
                   iterations = 5000,
                   replicates = 5,
                   temperature = 0.999,
                   diagnostic=FALSE)


par(mar=rep(5,4))
hist(unlist(lapply(resultsQuebec,function(x) return(x$score))))
qqline(unlist(lapply(output,function(x) return(x$score))))

writeResults<-function(output)
{


pvalues<-unlist(lapply(resultsQuebec,function(x) {stat<-x$pvalue;if(is.null(stat)) stat<-NA; return(stat)}))
subnetSize<-unlist(lapply(resultsQuebec,function(x) {stat<-x$size;if(is.null(stat)) stat<-NA; return(stat)}))
netSize<-unlist(lapply(resultsQuebec,function(x) {stat<-length(x$table$gene);if(is.null(stat)) stat<-NA; return(stat)}))

subnetScore<-unlist(lapply(resultsQuebec,function(x) {stat<-x$score;if(is.null(stat)) stat<-NA; return(stat)}))

abline(lm(subnetSize~netSize))



namesKegg<-names(pathways("hsapiens","kegg"))
}


iter<-1000
test<-array(NA,iter)
for (i in 1:iter)
{
  newscores<-data.frame(gene=scores$gene,score=sample(scores$score))
  sc<-searchSubnet(keggPathways[[sample(1:244,1)]],newscores,
                   quebecNull,
                   replicates=1,iterations=2000,
                   temperature=0.99,verbose=FALSE)$score
  print(i)
  if(!is.null(sc)) test[i]<-sc
}
hist(test,breaks=20)


library(qvalue)

subnetScore[
  which(qvalue(unlist(lapply(subnetScore,function(x) return(mean(test>x,na.rm=TRUE)))))$lfdr<0.05)
  ]



namesKegg[qvalue(unlist(lapply(subnetScore,function(x) return(mean(test>x,na.rm=TRUE)))),fdr.level=0.001)$significant]



namesKegg[which(qvalue(unlist(lapply(subnetScore,function(x) return(mean(test>x,na.rm=TRUE)))))$lfdr<0.05)]
subnetScore[which(qvalue(unlist(lapply(subnetScore,function(x) return(mean(test>x,na.rm=TRUE)))))$lfdr<0.05)]
subnetSize[which(qvalue(unlist(lapply(subnetScore,function(x) return(mean(test>x,na.rm=TRUE)))))$lfdr<0.05)]

which(qvalue(unlist(lapply(subnetScore,function(x) return(mean(test>x,na.rm=TRUE)))))$lfdr<0.05)[29]

(lapply(resultsQuebec,function(x) {stat<-x$table[x$table$state,]$gene;if(is.null(stat)) stat<-NA; return(paste(stat,sep=";"))}))

test[which(qvalue(unlist(lapply(subnetScore,function(x) return(mean(test>x,na.rm=TRUE)))))$lfdr<0.05)]


quebecKegg[[2]]
kegg[[2]]
cat(quebecKegg[[2]]$table[quebecKegg[[2]]$table$state,]$gene,sep="+")





### procedure overlapping
library(graphite)
?pathways
react<-pathways("hsapiens","reactome")
kegg<-pathways("hsapiens","kegg")
biocarta<-pathways("hsapiens","biocarta")
panther<-pathways("hsapiens","panther")
nci<-pathways("hsapiens","nci")


#pathway object to graphNEL
react<-convertIdentifiers(react, "entrez")
reactGraph<-lapply(react,pathwayGraph,edge.types=NULL)

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
save(overlap,file="overlapReactome.rda")

diag(overlap)<-0

sum(overlap>0.05,na.rm=TRUE)
mat<-matrix(NA,ncol=6,nrow=length(reactGraph)^2)
for(i in 1:length(reactGraph))
{
  for(j in 1:length(reactGraph))
  {
    line<-c(i,j,length(reactGraph[[i]]@nodes),length(reactGraph[[j]]@nodes),
            max(c(overlap[i,j],overlap[j,i])),
            sum(reactGraph[[i]]@nodes %in% reactGraph[[j]]@nodes))
    mat[i+(j-1)+(i-1)*length(reactGraph),]<-line
  }
  print(i)
}

mat2<-mat[mat[,5]>0.1&!is.na(mat[,5]),]

#case1
case1<-mat2[mat2[,5]==1,]
case1<-case1[case1[,2]%in%as.numeric(names(table(case1[,2])[table(case1[,2])>1])),]

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
