performK=function(X,warp_c, K, Ks,simu) 
{
N=dim(X)[2]
PCA=eigen(1/T*t(X)%*%(X))
#PCA=eigen(1/T*t(X)%*%(X))
PCAexp=sum(PCA$values[1:K])/sum(PCA$values)

#Adding scores of warping functions
PCA2=eigen(t(warp_c)%*%warp_c)
PCA2exp=sum(PCA2$values[1:Ks])/sum(PCA2$values)


if(Ks>0)
{	
	##After getting the "best" Amplitude space with K components we truncate the W space by Ks components
	cdec=eigen(t(warp_c)%*%warp_c)
	c_low=warp_c%*%cdec$vectors[,1:Ks]%*%t(cdec$vectors[,1:Ks])
	#print(dim(-c_low) )
	htinv=sapply(1:N, function(i)  warph(-c_low[,i]) ) 		###compute the inverse truncated warping function
	
	##Compute the K dimensional warping space
	eigDC=eigen(t(X)%*%(X))
	X_K=X%*%eigDC$vectors[,1:K]%*%t(eigDC$vectors[,1:K])  
	Xflow<<-sapply(1:N, function(i) approxfun(  t2,X_K[,i],rule = 2,method="linear"))
	X= sapply(1:N, function(i)  simu[[i]]$y) 
	Xint<<-sapply(1:N, function(i) approxfun(  t2,X[,i],rule = 2,method="linear"))
	Xt= sapply(1:N, function(i)  Xint[[i]](t2) )
	#print(dim(htinv) )
	Xflowat= sapply(1:N, function(i)  Xflow[[i]](htinv[,i]) )

	matplot(t,Xt, typ="l", col="black")
	#matlines(t,Xflowa, typ="l", col="grey")
	matlines(t,Xflowat, typ="l", col="grey")
	###No we "warp back" to compute the ResV= || X- X_K(h^-1_Ks(h(t)) ||
	residual =(1-mean(mean((Xflowat-Xt)^2))/mean(mean((X)^2)) )
	#ResVtruncated=
}
else
{
	eigDC=PCA
	X= sapply(1:N, function(i)  simu[[i]]$y) 
	X_K=X%*%eigDC$vectors[,1:K]%*%t(eigDC$vectors[,1:K])  
	residual=1-mean(mean((X-X_K)^2))/ mean(mean((X)^2))
	matplot(t,X, typ="l", col="black")
	#matlines(t,Xflowa, typ="l", col="grey")
	matlines(t,X_K, typ="l", col="grey")
	PCA2exp=0
}

list(Ampsp=PCAexp,Warpsp=PCA2exp, Explained=residual )
}
###To get a best model isterate over K and Ks (could be parallized) and choose the model where ResV is smallest


perform=function(X,c, kappa,P) 
{
	
kids=c(rep(0,39),rep(1,54))
PCA=eigen(1/T*t(X)%*%(X))
#PCA=eigen(1/T*t(X)%*%(X))
PCAexp=sum(PCA$values[1:kappa])/sum(PCA$values)

#Adding scores of warping functions
PCA2=eigen(t(c)%*%c)
PCA2exp=sum(PCA2$values[1:P])/sum(PCA2$values)

res=NULL
for(i in 1:93)
{
if(P==0)
{
	mylogit<-glm(kids[-i]~PCA$vectors[-i,1:kappa],family = "binomial")
	res[i]=1/(1+exp  (- c( 1, PCA$vectors[i,1:kappa]) %*% mylogit$coef) ) 
}
	
if(kappa==0)
{
	mylogit<-glm(kids[-i]~PCA2$vectors[-i,1:P],family = "binomial")
	res[i]=1/(1+exp  (- c( 1,PCA2$vectors[i,1:P]) %*% mylogit$coef) ) 
}	

if(P>0 && kappa>0)
{
	mylogit<-glm(kids[-i]~cbind( PCA$vectors[-i,1:kappa], PCA2$vectors[-i,1:P]),family = "binomial")
	res[i]=1/(1+exp  (- c( 1, Re(PCA$vectors[i,1:kappa]), PCA2$vectors[i,1:P]) %*% mylogit$coef) ) 
}	

}
false=1*(res<0.5)
true=1*(res>=0.5)
correct=sum(c(false[1:39],true[40:93]))
correct/93

if(P>0)
{	
	residual= performK(X,c,kappa,P)
}
else
{
	eigDC=eigen(t(X)%*%(X))
	X_K=X%*%eigDC$vectors[,1:kappa]%*%t(eigDC$vectors[,1:kappa])  
	residual=1-mean(mean((Xn-X_K)^2))/ mean(mean((Xn)^2))
	matplot(t,Xn, typ="l", col="black")
	#matlines(t,Xflowa, typ="l", col="grey")
	matlines(t,X_K, typ="l", col="grey")
}

list(MSE=mean( (kids-res)^2), Correct=correct/93,Ampsp=PCAexp,Warpsp=PCA2exp, Explained=residual[3])
}
