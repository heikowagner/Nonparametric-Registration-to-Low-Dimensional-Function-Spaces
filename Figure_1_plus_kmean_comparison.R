rm(list= ls())
library(fda)
library(trust)
library(fdasrvf)


#setwd("I:/OwnCloud2/Uni/Registration/R  code geordnet paper - Kopie")
#setwd("D:/Auslagerung/Dropbox/Dropbox/Registration/R code geordnet paper")
#setwd("D:/Dropbox/Dropbox/Registration/R code geordnet paper")
#setwd("I:/OwnCloud2/Uni/Programs/Registration")
#setwd("C:/Users/Heiko/Dropbox/Registration/R code geordnet paper")
#setwd("C:/Users/Heiko/Desktop/Registration")
setwd("D:/owncloud2/Reg_aufgeräumt")
setwd("I:/OwnCloud2/Reg_aufgeräumt")

set.seed(123)
#source("simulations/simu3.r", echo=TRUE, max.deparse.length=150)
source("simulations/kmean.r", echo=TRUE, max.deparse.length=150)
source("base/multireg2.r", echo=TRUE, max.deparse.length=150)


N=length(simu)
T=128
t=seq(from=0,to=1, length=T)
Xn=sapply(1:N, function(i) simu[[i]]$y)
Xn2=sapply(1:N, function(i) simu2[[i]]$y)

#Register curves using FR-Metric 
Freg=align_fPCA(Xn, t )

#Register curves using our Method
N=length(simu)
T1=length(simu[[1]]$y)
ts=seq(from=0,to=1,length=T1)
T=128
t=seq(from=0,to=1,length=T)
B=9

c02=rnorm(N*B,0,0.001)
outputK2=multireg2(simu,K=2,w2=0.0001,B=B, plot=0,c=c02, iter=20000, T2=T)

eigen_unreg= eigen( 1/dim(Xn)[1]* Xn%*%t(Xn))$values
eigen_K2= eigen( 1/dim(outputK2$X)[1]* outputK2$X%*%t(outputK2$X))$values
eigen_K1= eigen( 1/dim(Freg$fn)[1]* Freg$fn%*%t(Freg$fn))$values

##Plots for figure 2
par(mfrow=c(2,3))
matplot(t,Xn, typ="l" , xlab="t" , ylab="Y(t)", main="Unregisted model")
matplot(t,outputK2$X, typ="l" , xlab="t" , ylab="Y(h(t))", main="Registration with K=2")
matplot(t,Freg$fn, typ="l" , xlab="t" , ylab="Y(h(t))", main="FR-Metric Registration")

plot(log( eigen_unreg[1:5] ) , typ="b" , ylab="log(eigenvalue)", xlab="eigenvector", main="log-Eingenvalue comparison")
lines(log( eigen_K2[1:5] ) , typ="b", col="red")
lines(log( eigen_K1[1:5] ) , typ="b", col="blue")
legend(3.2,3, # places a legend at the appropriate place 
	c("No Registration","Registration with K=2","Registration using FR-Metric"), # puts text in the legend 
	lty=c(1,1,1), # gives the legend appropriate symbols (lines)
lwd=c(2.5,2.5,2.5),col=c("black","red","blue")) # gives the legend lines the correct color and width

matplot(t,outputK2$h, typ="l", main="Warpingfunctions", ylab="h", xlab="t")
matplot(t,Freg$gam, typ="l", main="Warpingfunctions", ylab="h", xlab="t")

###Comparison with K-Mean Cluster

#install.packages("fdakma")
library("fdakma")

###k-mean cluster
kmean=kma(x=t(t),y0=t(Xn), n.clust = 2 ,similarity.method="d0.pearson")
kma.show.results(kmean)

###New approach
par(mfrow=c(3,2))
matplot(t,Xn, typ="l" , xlab="t" , ylab="Y(t)", main="Unregisted model")
matplot(t,Xn2, typ="l" , xlab="t" , ylab="Y(h(t))", main="2 Dimensional true model")
#matplot(t,Freg$fn, typ="l" , xlab="t" , ylab="Y(t)", main="FR Metric Registration")
matplot(t,outputK2$X, typ="l" , xlab="t" , ylab="Y(t)", main="New Algorithm")
matplot(t(kmean$x.final),Xn, typ="l" , xlab="t" , ylab="Y(t)", main="K-means using fdakma Package")
matplot( svd(outputK2$X)$v[,1:2]  , main="First two scores a_ij")
plot(kmean$labels, main="Clusters identified by k-mean")




