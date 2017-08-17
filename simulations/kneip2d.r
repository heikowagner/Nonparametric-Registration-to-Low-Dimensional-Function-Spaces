#Combining Registration and Fitting for Functional Models
##########################################Simulating Data
set.seed(100)
n<-26
T=256
#ystar function
#z<-rnorm(2*n,1,0.25)

z1<-rnorm(n,1,0.35)
z2<-rnorm(n,1,0.6)


ystar<-function(h,z1,z2,i){


#z[2*i-1]*exp(-(h-1.5)^2/2) + z[2*i]*exp(-(h+1.5)^2/2)
z1[i]*exp(-(h-1.5)^2/2) + z2[i]*exp(-(h+1.5)^2/2)

}

#warping function
h<-function(t,a,i){


if (a[i]==0) {t}

else {6*(expm1(a[i]*(t+3)/6))/(expm1(a[i]))-3}


}

simu<-vector("list",n)

simu2<-vector("list",n)

hfkt<-vector("list",n)

t<-seq(-3,3,length.out=T)
#a<-seq(-1,1,length.out=n)
a<-rnorm(n)/2
dat<-NULL

for(i in 1:n){

simu[[i]]<-list(x=seq(-3,3,length.out=T),y=( ystar(h(t,a,i),z1,z2,i) ) )

#simu[[i]]<-list(x=seq(-3,3,length.out=101),y=ystar(h(t,a,i),z1,z2,i) )

class(simu)="curves"

simu2[[i]]<-list(x=seq(-3,3,length.out=T),y=ystar(t,z1,z2,i))
class(simu2)="curves"

#simu[[i]]=matrix(x=seq(-3,3,length.out=101),y=ystar(h(t,a,i),z,i))
#dat<-cbind(dat,ystar(h(t,a,i),z,i))
hfkt[[i]]=h(t,a,i)
}