rm(list= ls())


setwd("D:/owncloud2/Reg_aufgeräumt")
setwd("I:/OwnCloud2/Reg_aufgeräumt")
source("simulations/simu2c.r", echo=TRUE, max.deparse.length=150)
source("base/multireg2.r", echo=TRUE, max.deparse.length=150)

B=7
N=length(simu)
c02=rnorm(N*(B+1),0,0.001)

output8=multireg2(simu,K=8,w2=0.001,B=B, plot=0,c=c02, iter=15000)
output7=multireg2(simu,K=7,w2=0.001,B=B, plot=0,c=c02, iter=15000)
output6=multireg2(simu,K=6,w2=0.001,B=B, plot=0,c=c02, iter=15000)
output5=multireg2(simu,K=5,w2=0.001,B=B, plot=0,c=c02, iter=15000)
output4=multireg2(simu,K=4,w2=0.001,B=B, plot=0,c=c02, iter=15000)
output3=multireg2(simu,K=3,w2=0.001,B=B, plot=0,c=c02, iter=15000)
output2=multireg2(simu,K=2,w2=0.001,B=B, plot=0,c=c02, iter=15000)
output1=multireg2(simu,K=1,w2=0.001,B=B, plot=1,c=output3$c_pass, iter=15000)

#output1=multireg2(simu,K=1,w2=0.3,B=B, plot=1,c=output3$c,iter=5000,T2=64)

#output1a=multireg2(simu,K=1,w2=0.01,B=B, plot=1,c=c02,T2=64,iter=3000)

#test8= newuoa(scorew,getlam2,L=8,control = list(maxfun=3000))
#test7= newuoa(scorew,getlam2,L=7,control = list(maxfun=3000))
#test6= newuoa(scorew,getlam2,L=6,control = list(maxfun=3000))
#test5= newuoa(scorew,getlam2,L=5,control = list(maxfun=3000))
#test4= newuoa(scorew,getlam2,L=4,control = list(maxfun=3000))
#test3= newuoa(scorew,getlam2,L=3,control = list(maxfun=3000))
#test2= newuoa(scorew,getlam2,L=2,control = list(maxfun=3000))
#test1= newuoa(test3$par,getlam2,L=1,control = list(maxfun=9000))
##test32= newuoa(test3$par,getlam2,L=2,control = list(maxfun=3000))
#test1= newuoa(test3$par,getlam2,L=1,control = list(maxfun=3000))


K=c(29,8,7,6,5,4,3,2,1)
W=c( 0, var(c(output8$c)),var(c(output7$c)),var(c(output6$c)),var(c(output5$c)),var(c(output4$c)),var(c(output3$c)),var(c(output2$c)),var(c(output1$c)))

t2=seq(from=0,to=1, length=128)
t3=seq(from=0,to=1, length=64)
Xn=sapply(1:29,function(i) simu[[i]]$y)

par(mfrow=c(2,3))
matplot(t[106:406],Xn[106:406,], typ="l", ylab="X", xlab="t", main="No Registration (K=29)")
matplot(t2[26:101],output3$X[26:101,], ylab="X", xlab="t", typ="l", main="Registration to K=3 model")
matplot(t2[26:101],output1$X[26:101,], ylab="X", xlab="t", typ="l", main="Registration to K=1 model")
plot(K,W,typ="b" , ylab="var(W)" ,main="Var(W) depending on choice of K")
matplot(t2,output3$h, ylab="h", xlab="t",typ="l", main="Warpingfunctions for K=3 model")
matplot(t2,output1$h, ylab="h",  xlab="t", typ="l", main="Warpingfunctions for K=1 model")
t