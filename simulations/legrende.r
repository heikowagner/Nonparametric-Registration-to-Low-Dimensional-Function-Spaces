###Verschobene Normalverteilung
#install.packages("orthopolynom")
#library(orthopolynom)


#normalized.p.list <- legendre.polynomials( 10, normalized=TRUE )
#print( normalized.p.list )
set.seed(123)
N=15
a=  rnorm(N) 
b=  rnorm(N)  

t=seq(from=-1, to=1,by=0.02)

#f1= dnorm(t+0.8)%*%t(a) 
#f2= dnorm(t-0.8)%*%t(b)


f1=(-0.7905694 + 2.371708*t^2)%*%t(a)  
#f1=(4.397265*t - 20.52057*t^3 + 18.46851*t^5)%*%t(a) 

f2=(0.7954951 - 7.954951*t^2 + 9.280777*t^4)%*%t(b) 

#f2=(0.7954951 - 7.954951*t^2 + 9.280777*t^4)%*%t(b)


f=-( f1+f2 )

matplot(f, typ="l")



h<-function(t,a,i){


if (a[i]==0) {t}

else {
	#6*(expm1(a[i]*(t+3)/6))/(expm1(a[i]))-3
	2*(expm1(a[i]*(t+1)/2))/(expm1(a[i]))-1
	
	}


}


a<-seq(-1,1,length.out=N)
simu=NULL
simuorig=NULL
for(i in 1:(N)){

#simu[[i]]<-list(x=t ,y=(approx(h(t,a,i),f[,i],t)$y+rnorm(t)*0.1+ rnorm(1)) )

simu[[i]]<-list(x=t ,y=(approx(h(t,a,i),f[,i],t)$y) )
class(simu)="curves"

simuorig[[i]]<-list(x=t ,y=f[,i] )
class(simu)="curves"

}
