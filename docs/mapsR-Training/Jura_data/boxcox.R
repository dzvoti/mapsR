

BoxCox<-function(x,lims){
library(MASS)
# estimate Box-Cox transformation parameters for values in x

if(missing(lims)){lambda<-seq(-5,5,0.001)}else{lambda<-seq(lims[1],lims[2],0.001)}

bc<-boxcox(x~1,lambda=lambda,plotit=F)

# Extract ML estimate of parameter from output

maxlik<-max(bc$y)
lam<-bc$x[which(bc$y==maxlik)]

# Find bounds of 95% confidence interval of BC parameter

CI95<-(maxlik-qchisq(0.95,1))
LB95<-bc$x[min(which(bc$y>CI95))]
UB95<-bc$x[max(which(bc$y>CI95))]



# Apply transformation to data

if(lam==0){z<-log(x)}else{
z<-(x^lam-1)/lam}


#Output plots

par(mfrow=c(2,2))
qqnorm(x,main="Raw data",pch=16)
qqline(x,lwd=2,col="red")
qqnorm(z,main="Transformed data",pch=16)
qqline(z,lwd=2,col="red")
hist(x,main="Raw data",col="AliceBlue")
hist(z,main="Transformed data",col="AliceBlue")

dev.new()

plot(bc$x,bc$y,type="l",xlab=expression(lambda),ylab="Likelihood")

lines(c(lam,lam),c(min(bc$y),max(bc$y)),lty=5)
lines(c(UB95,UB95),c(min(bc$y),max(bc$y)),lty=3)
lines(c(LB95,LB95),c(min(bc$y),max(bc$y)),lty=3)

return(z)
}