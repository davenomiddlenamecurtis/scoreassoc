# setwd("/Users/dave_000/OneDrive/sharedseq/SSS2/scores/SSS2.vrare.1/combinedscores")
setwd("d:/dave/OneDrive/sharedseq/SSS2/scores/SSS2.vrare.1/combinedscores")
scores.0=read.table("scores.fifth.0.out",header=FALSE)
scores.1=read.table("scores.fifth.1.out",header=FALSE)
scores.2=read.table("scores.fifth.2.out",header=FALSE)
scores.3=read.table("scores.fifth.3.out",header=FALSE)
scores.4=read.table("scores.fifth.4.out",header=FALSE)

scores=merge(scores.0,scores.1,all=TRUE)
scores=merge(scores,scores.2,all=TRUE)
scores=merge(scores,scores.3,all=TRUE)
scores=merge(scores,scores.4,all=TRUE)
attach(scores)
t.test(V4[V2==1],V4[V2==0],var.equal=TRUE)
scores$cc<-with(scores,ifelse(V2==1,as.character("CASE"),as.character("CONTROL")))
# because having numeric value does not work
library(ggplot2)
ggplot(scores,aes(x=V4,colour=cc))+
	geom_density()
library(pROC)
roc(scores$V2,scores$V4,levels=c(0,1))

scores.minimal.0=read.table("scores.minimal.fifth.0.out",header=FALSE)
scores.minimal.1=read.table("scores.minimal.fifth.1.out",header=FALSE)
scores.minimal.2=read.table("scores.minimal.fifth.2.out",header=FALSE)
scores.minimal.3=read.table("scores.minimal.fifth.3.out",header=FALSE)
scores.minimal.4=read.table("scores.minimal.fifth.4.out",header=FALSE)
scores.minimal=merge(scores.minimal.0,scores.minimal.1,all=TRUE)
scores.minimal=merge(scores.minimal,scores.minimal.2,all=TRUE)
scores.minimal=merge(scores.minimal,scores.minimal.3,all=TRUE)
scores.minimal=merge(scores.minimal,scores.minimal.4,all=TRUE)
attach(scores.minimal)
t.test(V4[V2==1],V4[V2==0],var.equal=TRUE)
scores.minimal$cc<-with(scores.minimal,ifelse(V2==1,as.character("CASE"),as.character("CONTROL")))
# because having numeric value does not work
library(ggplot2)
ggplot(scores.minimal,aes(x=V4,colour=cc))+
	geom_density()
roc(scores.minimal$V2,scores.minimal$V4,levels=c(0,1))

