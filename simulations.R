### Example
require(graphite);require(graph);require(igraph)


i<-2
ii<-2
iii<-2
iiii<-3
iiiii<-1

  netSize<-c(20,50,100)[i]
  subnetSize<-c(5,10,20)[ii]
  connectivity<-c(0.2,0.5,1)[iii]
  highScores<-c(2,3,5)[iiii]
  iter<-c(5000,10000,50000)[iiiii]
  cComp<-5

  filename<-paste(i,ii,iii,iiii,iiiii,iiiiii,"output.rda",sep="_")

  OK<-FALSE
  while(!OK) {
  g1 <- randomEGraph(paste(1:netSize),0.7)
  Len<-unlist(lapply(graph::connComp(g1),length))
  if(max(Len)<subnetSize){
    OK<-FALSE
  } else {
    OK<-TRUE
  }
  }
  X<-unlist(graph::connComp(g1)[(Len==max(Len))[1]])

  connected_comp<-graph::subGraph(X,g1)
  g1<-graph::ugraph(connected_comp)

  scores<-data.frame(genes=nodes(g1),scores=(rnorm(length(nodes(g1)))))
  bkgd<-data.frame(k=1:netSize,mu=rep(0,netSize),sigma=rep(1,netSize))

  trueVal<-c(1,rep(NA,subnetSize-1))
  for(i in 2:subnetSize){
    adjacent<-unlist(adj(g1,as.character(trueVal)))
    adjacent<-unique(adjacent[!adjacent%in%trueVal])
    trueVal[i]<-sample(adjacent,1)
  }

  scores[trueVal,2] <- rnorm(subnetSize,mean=highScores)

  Rprof(filename = "Rprof.out",line.profiling=TRUE)
  signetObject <- searchSubnet(pathway = g1,
                               nullDist = bkgd,
                               scores = scores,
                               subnetScore = "ideker",
                               iterations = 1000,
                               replicates=1)
  Rprof(NULL)
  summaryRprof(filename = "Rprof.out")

  TP<-sum(signetObject$network$gene[signetObject$network$active]%in%trueVal)
  FP<-sum(!signetObject$network$gene[signetObject$network$active]%in%trueVal)
  TN<-sum(!signetObject$network$gene[!signetObject$network$active]%in%trueVal)
  FN<-sum(signetObject$network$gene[!signetObject$network$active]%in%trueVal)

  output<-c(netSize=netSize,
            subnetSize=subnetSize,
            iter=iter,
            connectivity=connectivity,
            highScores=highScores,
            replicate=replicate,
            TP=TP,FP=FP,TN=TN,FN=FN)

  save(signetObject,output,file=filename)

par(mar=c(3,3,3,3))
plot(OUT$FP/OUT$netSize,OUT$TP/OUT$netSize)

# get from cluster
OUT<-NULL
for(i in 1:3){
  for(ii in 1:3){
    for(iii in 1:3){
      for(iiii in 1:3){
        for(iiiii in 1:3){
          for(iiiiii in 1:20){
            filename<-paste("//cmpgmatrix.unibe.ch/gouy/simulations/archive/",i,"_",ii,"_",iii,"_",iiii,"_",iiiii,"_",iiiiii,"_output.rda",sep="")
            try(load(filename))
            OUT<-rbind(OUT,output)
          }
        }
      }
    }
  }
}

#with mu0 qnd 1
for(i in 4:5){
  for(ii in 1:200){
    filename<-paste("//cmpgmatrix.unibe.ch/gouy/simulations/archive/3_2_2_",i,"_1_",ii,"_output.rda",sep="")
    try(load(filename))
    OUT<-rbind(OUT,output)  }
}

rownames(OUT)<-NULL
OUT<-as.data.frame(OUT)
OUT$FPR<-OUT$FP/(OUT$FP+OUT$TN)
OUT$PPV<-OUT$TP/(OUT$TP+OUT$FP)
OUT$TPR<-OUT$TP/(OUT$TP+OUT$FN)
OUT$F1<-2*OUT$TP/(2*OUT$TP+OUT$FN+OUT$FP)
OUT$correct<-OUT$FP
OUT$correct[OUT$correct==0]<-0
OUT$correct[!OUT$correct==0]<-1

# save(OUT,file="simulations.rda")
load(file="simulations.rda")


# boxplot(OUT$FPR~(OUT$highScores),ylim=c(0,1))
m<-lm(OUT$FP~OUT$netSize+OUT$subnetSize+OUT$highScores+OUT$connectivity+OUT$iter)
boxplot(OUT$FPR~OUT$subnetSize)
boxplot(OUT$FPR~OUT$highScores)

m<-glm(OUT$correct~OUT$netSize+OUT$subnetSize+OUT$highScores+OUT$connectivity+OUT$iter,family = "binomial")

m<-glm(cbind((OUT$TP),(OUT$FP))~OUT$netSize+OUT$subnetSize+OUT$highScores+OUT$connectivity+OUT$iter,family = "binomial")
summary(m)
anova(m)
# m<-glm(cbind((OUT$TP+OUT$TN),(OUT$FP+OUT$FN))~OUT$netSize*OUT$subnetSize*OUT$highScores*OUT$connectivity*OUT$iter,family = "binomial")
# summary(m)
# anova(m)

require(blackbox)
m<-glm(cbind((TP+TN),(FP+FN))~highScores+netSize+subnetSize,family = "binomial",data=OUT)
summary(m)
anova(m)

m<-glm(correct~highScores+netSize+subnetSize,family = "binomial",data=OUT)
summary(m)
anova(m)

grid<-init_grid(lower=c(highScores = 0, netSize = 20, subnetSize = 5),
          upper=c(highScores = 6, netSize = 100, subnetSize = 30),
          nUnique = 500)

yweight <- predict(m, grid,type="response")

par(mfrow=c(1,3))

plot(grid$highScores, 1-yweight,pch=16,col=rgb(0,0,0,alpha=0.3),cex=0.8)
abline(h=c(0.1,0.05,0.01),lty=3)
plot(grid$netSize, 1-yweight,pch=16,col=rgb(0,0,0,alpha=0.3),cex=0.8)
abline(h=c(0.1,0.05,0.01),lty=3)
plot(grid$subnetSize, 1-yweight,pch=16,col=rgb(0,0,0,alpha=0.3),cex=0.8)
abline(h=c(0.1,0.05,0.01),lty=3)

#significance

head(OUT)
scor<-NULL
for(i in 1:3){
  for(ii in 1:3){
    for(iii in 1:3){
      for(iiii in 1:3){
        for(iiiii in 1:3){
          for(iiiiii in 1:20){
            filename<-paste("//cmpgmatrix.unibe.ch/gouy/simulations/archive/",i,"_",ii,"_",iii,"_",iiii,"_",iiiii,"_",iiiiii,"_output.rda",sep="")
            try(load(filename))
            scor<-rbind(scor,mean(signetObject$network$score[signetObject$network$active]))
          }
        }
      }
    }
  }
}

plot(dnorm(scor),OUT$FP)
pval <- 1-pnorm(scor)
hist(pval)

plot(pval,pval2)

abline(0,1)
abline(h=0.05,v=0.05,col=2,lty=2)
plot(scor,pnorm(scor, mean = 0, sd = 1,lower.tail = FALSE))

# new pvalues with empirical distribution
scor2<-NULL
for(i in 1:1500){
  filename<-paste("//cmpgmatrix.unibe.ch/gouy/simulations/null_dist/",i,"_null.rda",sep="")
  (load(filename))
  scor2<-rbind(scor2,mean(signetObject$network$score[signetObject$network$active]))
}
hist(scor2)

pval2 <- unlist(lapply(scor,function(x)mean(x<scor2)))
hist(pval2)
?cut
pval2[pval2==0]<-1e-3
pclass<-cut(log10(pval),breaks=0:-10)

tapply(OUT$FPR[cond],pclass[cond],mean)
pclass2<-cut(log10(pval2),breaks=0:-10)
tapply(OUT$FPR[cond],pclass2[cond],mean)

cond<-which(OUT$highScores>0)

## fixed parameters
subs <- OUT[which(OUT$highScores==2),]
hist(subs$FPR)
subs <- OUT[which(OUT$highScores==3),]
hist(subs$FPR)
subs <- OUT[which(OUT$highScores==5),]
hist(subs$FPR)

FP

subs <- OUT[which(OUT$netSize == 100 & OUT$subnetSize == 10),]
subs <- OUT[which(OUT$highScores == 5 & OUT$subnetSize == 10),]
subs <- OUT[which(OUT$highScores == 5 & OUT$netSize == 100),]

hist(subs$FPR)
boxplot(subs$FP~subs$subnetSize, xlab="subnet size",ylab="FP")
abline(h=0.05)

subs <- which(OUT$netSize == 100 & OUT$subnetSize == 10)
boxplot(pval2[subs]~OUT$highScores[subs])
subs <- which(OUT$highScores == 5 & OUT$subnetSize == 10)
boxplot(pval2[subs]~OUT$netSize[subs])
subs <- which(OUT$highScores == 5 & OUT$netSize == 100)
boxplot(pval2[subs]~OUT$subnetSize[subs])

boxplot(pval2~OUT$highScores,ylim=c(0,0.25))

hist(pval*4)
hist(pval2)

subs <- which(OUT$netSize == 100 & OUT$subnetSize == 10 & OUT$highScores == 3)
plot(OUT$FP[subs],pval2[subs],main="mu=3")

subs <- which(OUT$netSize == 100 & OUT$subnetSize == 10 & OUT$highScores == 2)
plot(OUT$FP[subs],pval2[subs],main="mu=2")

subs <- which(OUT$netSize == 100 & OUT$subnetSize == 10 & OUT$highScores == 5)
plot(OUT$FP[subs],pval2[subs],main="mu=5")




#@plot a run
par(mfrow=c(1,1))
matplot(signetObject$simulated_annealing[,-2],cex=0.1,type="l",lty=1)





pdf(height = 10,width = 2,"p-values distributions.pdf")
par(mfrow=c(5,1))
for(i in c(0,1,2,3,5)){
  subs<-which(OUT$highScores==i)
  hist(pval[subs],breaks=20,xlim=c(0,0.3))
}
for(i in c(0,1,2,3,5)){
  subs<-which(OUT$highScores==i)
  hist(pval2[subs],breaks=20,xlim=c(0,1))
}
dev.off()

par(mfrow=c(1,1))
boxplot((OUT$FP+OUT$FN)/OUT$netSize ~ OUT$highScores)


mean(res[[1]]$table[res[[1]]$table$state,]$score)
hist(unlist(lapply(res,function(x) mean(x$table[x$table$state,]$score))))


#TP/FP+TP
tapply((OUT$FPR),cut(log10(pval2),breaks=0:-4),mean)

#pval2 FPR
plot((OUT$FP)/OUT$netSize ~ log10(pval2))
#score FPR
par(mfrow=c(1,2))
boxplot((OUT$FP)/OUT$netSize ~ OUT$highScores,ylim=c(0,0.5))
subs<-which(pval2<0.01)
boxplot((OUT$FP[subs])/OUT$netSize[subs] ~ OUT$highScores[subs],ylim=c(0,0.5))


par(mfrow=c(1,2))
boxplot((OUT$TP)/(OUT$FP+OUT$TP) ~ OUT$highScores,ylim=c(0,1))
subs<-which(pval2<0.01)
boxplot((OUT$TP[subs])/(OUT$FP[subs]+OUT$TP[subs]) ~ OUT$highScores[subs],ylim=c(0,1))

(1-tapply((OUT$TP)/(OUT$FP+OUT$TP),OUT$highScores,mean))/0.9
#score without non significant FPR





#overlapping coefficient between two normal distributions as a function of the mean
OUT2<-OUT[which(OUT$netSize == 100 & OUT$subnetSize  == 10& OUT$iter == 5000),]
OUT2<-OUT
(OUT2$FP)/(2*OUT2$FP+2*OUT2$TP)/(1-OUT2$subnetSize/OUT2$netSize)
m<-glm(cbind((FP),(TP))~highScores,family = "binomial",
       data=OUT)


m<-glm(cbind((FP),(FP+TN+TP))~highScores+netSize+subnetSize,family = "binomial",
       data=OUT)



grid<-init_grid(lower=c(highScores = 0),
                upper=c(highScores = 6),
                nUnique = 500)

xa<-seq(0,6,0.1)
yweight <- predict(m, data.frame(highScores=xa,netSize=100,subnetSize=20),type="response")
y<-2*pnorm(-abs(xa)/2)
plot(xa, (yweight),pch=16,col=rgb(0,0,0,alpha=1),cex=0.8,
     ylim=c(0,0.5))
lines(xa,y/2,type="l",lwd=1,lty=2)
plot((OUT2$FP)/(2*OUT2$FP+2*OUT2$TP) ~ OUT2$highScores,pch=16,col=rgb(0,0,0,alpha=0.05))

ms<-tapply((OUT2$FP+OUT2$FN)/(OUT2$netSize), OUT2$highScores,mean)
sds<-tapply((OUT2$FP+OUT2$FN)/(OUT2$netSize), OUT2$highScores,sd)

points(c(0,1,2,3,5),ms,cex=2,pch=16,col=2)


points(OUT$highScores,OUT$PPV,pch=16,col=rgb(0,0,0,alpha=0.05))


#P(FP)
OUT2<-OUT[which(OUT$netSize == 100 & OUT$subnetSize  == 10),]
OUT2<-OUT
xa<-seq(0,6,0.1)
YOBS<-(OUT2$FP)/(2*OUT2$FP+2*OUT2$TP)/(1-OUT2$subnetSize/OUT2$netSize)
# 
# YOBS<-(OUT2$FP)/(OUT2$netSize)
# /(1-OUT2$subnetSize/OUT2$netSize)

plot(YOBS ~ OUT2$highScores,pch=16,col=rgb(0,0,0,alpha=0.05))


y2<-(2*pnorm(-abs(xa)/2))/(1+2*pnorm(-abs(xa)/2))
y3<-(pnorm(-abs(xa)/2))/(1+pnorm(-abs(xa)/2))

ms<-tapply(YOBS, OUT2$highScores,mean)
sds<-tapply(YOBS, OUT2$highScores,sd)

plot(c(0,1,2,3,5),ms,cex=2,pch=16,col=2,ylim=c(0,0.55))
segments(x0=c(0,1,2,3,5),y0=ms+sds,y1=ms-sds,col=2,lwd=2)
lines(xa,y2,lty=2,lwd=2)
m

Estimate Std. Error z value Pr(>|z|)    
(Intercept) -1.0179053  0.0274608  -37.07   <2e-16 ***
  highScores  -0.4638224  0.0042835 -108.28   <2e-16 ***

  
lines(xa,1-1/(1+exp(-1.01+-0.463822*xa)),col=2)
xa
