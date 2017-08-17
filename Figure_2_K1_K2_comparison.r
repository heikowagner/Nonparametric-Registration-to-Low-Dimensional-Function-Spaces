rm(list= ls())


#setwd("I:/OwnCloud2/Uni/Registration/R  code geordnet paper - Kopie")
#setwd("D:/Auslagerung/Dropbox/Dropbox/Registration/R code geordnet paper")
#setwd("D:/Dropbox/Dropbox/Registration/R code geordnet paper")
#setwd("I:/OwnCloud2/Uni/Programs/Registration")
#setwd("C:/Users/Heiko/Dropbox/Registration/R code geordnet paper")
#setwd("C:/Users/Heiko/Desktop/Registration")
setwd("D:/owncloud2/Reg_aufgeräumt")
setwd("I:/OwnCloud2/Reg_aufgeräumt")

source("base/multireg2.r", echo=TRUE, max.deparse.length=150)
#source("simulations/kneip2d.r", echo=TRUE, max.deparse.length=150)
source("simulations/legrende.r", echo=TRUE, max.deparse.length=150)

B=7
N=length(simu)
c02=rnorm(N*(B),0,0.001)


outputK2=multireg2(simu,K=2,w2=0.0001,B=B, plot=0,c=c02, iter=30000, T2=512)
outputK1=multireg2(simu,K=1,w2=0,B=B, plot=0,c=c02, iter=30000, T2=512)

Xn=sapply(1:N, function(i) simu[[i]]$y)

eigen_unreg= eigen( 1/dim(Xn)[1]* Xn%*%t(Xn))$values
eigen_K2= eigen( 1/dim(outputK2$X)[1]* outputK2$X%*%t(outputK2$X))$values
eigen_K1= eigen( 1/dim(outputK1$X)[1]* outputK1$X%*%t(outputK1$X))$values

par(mfrow=c(2,3))
matplot(t,Xn, typ="l", ylab="X", xlab="t", main="No Registration (K=29)")
matplot(t2,outputK2$X, ylab="X", xlab="t", typ="l", main="Registration to K=2 model")
matplot(t2,outputK1$X, ylab="X", xlab="t", typ="l", main="Registration to K=1 model")

plot(log( eigen_unreg[1:5] ) , typ="b", ylim=c(-3.5,2.5) , ylab="log(eigenvalue)", xlab="eigenvector", main="log-Eingenvalue comparison")
lines(log( eigen_K2[1:5] ) , typ="b", col="red")
lines(log( eigen_K1[1:5] ) , typ="b", col="blue")

legend(3.5,2, # places a legend at the appropriate place 
	c("No Registration","Registration with K=2","Registration with K=1"), # puts text in the legend 
	lty=c(1,1,1), # gives the legend appropriate symbols (lines)
lwd=c(2.5,2.5,2.5),col=c("black","red","blue")) # gives the legend lines the correct color and width


matplot(t2,outputK2$h, ylab="h", xlab="t",typ="l", main="Warpingfunctions for K=2 model")
matplot(t2,outputK1$h, ylab="h",  xlab="t", typ="l", main="Warpingfunctions for K=1 model")