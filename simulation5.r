###Registration only over the functional space via optimize
rm(list= ls())


library(fda)
library(trust)

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

N=length(simu)
T1=length(simu[[1]]$y)
ts=seq(from=0,to=1,length=T1)
T=512
t=seq(from=0,to=1,length=T)
B=7

c02=rnorm(N*(B+1),0,0.001)

#output4=multireg2(simu,K=2,w2=0.001,B=B, plot=1,c=c02, iter=15000)
#output3=multireg2(simu,K=3,w2=0.001,B=B, plot=0,c=c02, iter=5000)
output2=multireg2(simu,K=3,w2=0.0001,B=B, plot=1,c=c02, iter=9000, T2=T)
#output1=multireg2(simu,K=1,w2=0.01,B=B, plot=0,c=c02, iter=5000)
#test3= newuoa(scorew,getlam2,L=3,Xf=Xf3, control = list(maxfun=3000))


#Xn=sapply(1:N, function(i) simu[[i]]$y )
Xn=datx2[,sele]
###Determine Phase and Amplitude
tt=seq(0,1,length=T)

hk5=output2$h
Xreg5=output2$X

K=c(29,8,7,6,5,4,3,2,1)
W=c( 0, NB(test8$par),NB(test7$par),NB(test6$par),NB(test5$par),NB(test4$par),NB(test3$par),NB(test2$par),NB(test1$par))

par(mfrow=c(2,3))
matplot(ts[106:406],datx[106:406,sele], typ="l", main="No Registration " , xlab="t", ylab="x(t)")
matplot(ts[106:406],vectors[100:400,1:3]%*%a2[1:3,sele], typ="l", main="3 Dimensional decomposition" , xlab="t", ylab="x(t)")
matplot(t[106:406],Xreg5[106:406,], typ="l", main="Registration to K=3 model", xlab="t", ylab="x(h(t))")



matplot(ts[106:406],Xn[106:406,], typ="l", main="No Registration ", xlab="t", ylab="x(t)")
matplot(ts[106:406],vectors[106:406,1:3]%*%a[1:3,sele], typ="l", main="3 Dimensional decomposition", xlab="t", ylab="x(t)")
#matplot(t,Xreg1, typ="l", main="Registration to K=3 model", xlab="t", ylab="x(t)",lty=1)
matplot(ts[106:406],vectors[106:406,1:3], typ="l", main="First three eigenfunctions", xlab="t", ylab="",lty=1)





