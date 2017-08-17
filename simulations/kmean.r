#Combining Registration and Fitting for Functional Models
##########################################Simulating Data


###Sangalli 2010



#set.seed(128)
n<-48
T=128
#ystar function
#z<-rnorm(2*n,1,0.25)

#z1<-rnorm(n,1,0.6)
#z2<-rnorm(n,0.2,0.1)

#z1<-abs(rnorm(n,4,0.6))
#z2<-abs(rnorm(n,0.1,0.15))

#z1<-rnorm(n,0,0.8)
#z2<-rnorm(n,0,0.15)

z1<-rnorm(n,0,0.2)
z2<-rnorm(n,0,0.1)
z3<-rnorm(n,0,0.3)
z3[1]=0

t=seq(0,2*pi,length=T)
ts=seq(0,1,length=T)



#curves=(1+z1)*sin(z3+(1+z4)*t)+(1+z2)*sin((z3+(1+z4)*t)^2/(2*pi) )

simu=NULL
simu2=NULL

h<-function(i){
if (z3[i]==0) {ts}
else {
	#t-0.02*(i-1)
	#(expm1(z3[i]*(ts)))/(expm1(z3[i]))
	#(expm1(z3[i]*(ts)))/(expm1(z3[i]))
	
	(exp(z3[i]*(ts))-1)/(exp(z3[i])-1)
	#exp( (ts -1)*z3[i])
	#2*pi*(expm1(z3[i]*(t)/(2*pi)))/(expm1(z3[i]))
	}
}



for(i in 1:(n/2))
{
simu[[i]]<-list(x=t,y=(1.4+z1[i])*sin( 2*pi*h(i) )+(1+z2[i])*sin(2*pi*(h(i))^2 ))
simu2[[i]]<-list(x=t,y=(1.4+z1[i])*sin( 2*pi*h(1) )+(1+z2[i])*sin(2*pi*(h(1))^2 ) )
#simu2[[i]]<-list(x=t,y=(1+z1[i])*sin(h(1))+(1+z2[i])*sin(h(1)^2/(2*pi) )*1.05 )
}

for(i in (n/2+1):n)
{
#simu[[i]]<-list(x=t,y=(2.2+z1[i])*sin(h(i))+(-1.1+z2[i])*sin(h(i)^2/(2*pi) ) )
simu[[i]]<-list(x=t,y=(2.1+z1[i])*sin(2*pi*h(i))+(-1+z2[i])*sin(2*pi*(h(i))^2 ) )
simu2[[i]]<-list(x=t,y=(2.1+z1[i])*sin(2*pi*h(1))+(-1+z2[i])*sin(2*pi*(h(1))^2 ) )
#simu[[i]]<-list(x=t,y=(2.2+z1[i])*sin(h(i))+(-1.2+z2[i])*sin(h(i)^2/(2*pi) ) )
#simu2[[i]]<-list(x=t,y=(2.2+z1[i])*sin(h(1))+(-1.2+z2[i])*sin(h(1)^2/(2*pi) ) )
#simu2[[i]]<-list(x=t,y=(2.2+z1[i])*sin(h(1))+(-1.1+z2[i])*sin(h(1)^2/(2*pi) ) )
}
par(mfrow=c(1,2))

Xn=sapply(1:n,function(i) simu2[[i]]$y)
matplot(ts,Xn, typ="l")


Xn2=sapply(1:n,function(i) simu[[i]]$y)
matplot(ts,Xn2, typ="l")

hi=sapply(1:n,function(i) h(i))
matplot(ts,hi, typ="l")