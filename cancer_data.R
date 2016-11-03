source("https://bioconductor.org/biocLite.R")
biocLite("hgu133plus2.db")
library(SPIA)
data(colorectalcancer)

hist(top$P.Value)

zscores<-(1-qnorm(top$P.Value))
z<-(zscores-mean(zscores))/sd(zscores)
qqnorm(zscores)
qqline(zscores)


cancer.data<-data.frame(gene=top$ENTREZ,score=zscores)
head(cancer.data)
hist(cancer.data$score,breaks=25)

library(signet)
