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

quebecNull<-nullDistribution(keggPathways,scores = scores,kmin=2,kmax=50,iterations=100)
save(quebecNull,file="quebecNull.rda")
load("quebecNull.rda")
#before generating the entire null distribution
for(i in 1:450) quebecNull<-rbind(quebecNull,quebecNull[49,])
quebecNull$k<-2:500
plot(quebecNull)

quebecKegg<-lapply(keggPathways[1:50],searchSubnet,
                         scores=scores,
                         nullDist = quebecNull,
                         iterations = 5000,
                         replicates = 5,
                         temperature = 0.999,
                         diagnostic=FALSE)
quebecKegg2<-lapply(keggPathways[51:100],searchSubnet,
                   scores=scores,
                   nullDist = quebecNull,
                   iterations = 5000,
                   replicates = 5,
                   temperature = 0.999,
                   diagnostic=FALSE)
quebecKegg3<-lapply(keggPathways[100:150],searchSubnet,
                   scores=scores,
                   nullDist = quebecNull,
                   iterations = 5000,
                   replicates = 5,
                   temperature = 0.999,
                   diagnostic=FALSE)
quebecKegg4<-lapply(keggPathways[151:200],searchSubnet,
                    scores=scores,
                    nullDist = quebecNull,
                    iterations = 5000,
                    replicates = 5,
                    temperature = 0.999,
                    diagnostic=FALSE)
quebecKegg5<-lapply(keggPathways[201:244],searchSubnet,
                    scores=scores,
                    nullDist = quebecNull,
                    iterations = 5000,
                    replicates = 5,
                    temperature = 0.999,
                    diagnostic=FALSE)



quebecKegg[[2]]
kegg[[2]]
cat(quebecKegg[[2]]$table[quebecKegg[[2]]$table$state,]$gene,sep="+")
