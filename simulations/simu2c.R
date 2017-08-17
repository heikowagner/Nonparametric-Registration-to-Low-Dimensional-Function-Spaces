nobs<-29
shift<-seq(0.25,0.75,length.out=nobs)
sscale<-5*(shift-0.5)^2/2+1
#t<-seq(-0.9,0.9,length.out=512)
t<-seq(0,1,length.out=512)
daten<-vector("list", nobs)
datx<-NULL
orient=sample(c(-1,1), replace=TRUE,nobs)
orig<-vector("list", nobs)
for(j in 1:nobs){
#datx<-cbind(datx, dnorm(t,shift[j],0.05)*sscale[j])
#datx<-cbind(datx, orient[j]*dbeta((t-0.5)*1.5+shift[j],50,50)*sscale[j])
datx<-cbind(datx, dbeta((t-0.5)*1.5+shift[j],50,50)*sscale[j])
}
matplot(t,datx,type="l")

simu=NULL
for(i in 1:nobs){

simu[[i]]<-list(x=t ,y= datx[,i] )
class(simu)="curves"

}

#t<-seq(0,1,length.out=512)
#plot(dbeta(t,12,12))