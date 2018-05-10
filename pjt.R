set.seed(1234)
R = 2000 # Number of Gibbs iterations
tenit = R/10 # To be used later for the last 10% of the Gibbs iterations

library(ggplot2)
test1 <- read.csv("~/Documents/OneDrive/MATH640/test1.csv", sep="")
Y = test1$DAILY_AQI_VALUE
plot(density(Y))
empty_k = 0
Y = scale(Y)
N = length(Y)
plot(density(Y))
ggplot()+
  geom_density(aes(Y),fill = "green", alpha = 0.2)+
  xlim(c(-3,3))+
  ggtitle('The density of PM2.5 for Hawaii and LA')

summary(Y)

t = seq(0.1,20,0.1)
#% Now ready to perform Gibbs sampling
#% Prior variance of mu
sigma0 = .15
#% Prior Dirichlet concentration of indicator probability p
c = 1e-3# % 1e-4 for Test 12
#% Likelihood variance of mu
sigma = 0.3
#% Initialization of the variables discussed
Kr = rep(1,R)
Kr2 = rep(1,R)
xr = matrix(rep(1,R*N),nrow=R)
mur = matrix(nrow=R,ncol=N)
nr = matrix(nrow=R,ncol=N)
mur[1,1] = 0
for(r in 1:(R-1)){
  x = xr[r,]
  K = Kr[r]
  mu = mur[r,1:K]
  for(n in 1:N){
    pxy = rep(1,K+1)
    for(k in 1:K){
      i = (x==k)
      Nkni = sum(i)
      if(x[n]==k){
        Nkni = Nkni - 1
      }
    pxy[k] = Nkni*dnorm(Y[n],mu[k],sigma)
    }
  pxy[K+1] = c*dnorm(Y[n],0,sigma0+sigma)
  pxy = pxy/sum(pxy)
  
  x[n] = sample(K+1,1,prob = pxy)
  
  if(x[n]>K){
    K=K+1
    muhat = sigma0*Y[n]/(sigma+sigma0)
    sigmahat = sigma*sigma0/(sigma+sigma0)
    mu[K]=sqrt(sigmahat)*rnorm(1)+muhat
    }
  }
  
  
  Kr2[r] = K
  for(k in 1:K){
    i = (x==k)
    nr[r,k] = sum(i)
    if(sum(i)==0){
      Kr2[r] = Kr[r]-1
      mu[k] = 0
    }
    else{
    
    
    muhat = sigma0*sum(Y[i])/(sigma+sigma0*sum(i))
    sigmahat = sigma*sigma0/(sigma+sigma0*sum(i))
    mu[k] = sqrt(sigmahat)*rnorm(1)+muhat}
  }
  
  
  xr[r+1,]=x
  Kr[r+1]= K
  mur[r+1,1:K] = mu
  
  xrp = xr[(R-tenit):R,]
  murp = mur[(R-tenit):R,]
  
  xmode = mode(xrp)
  mumean = mean(murp)
}




x = 1:2000
dmur=cbind(x,mur)
colnames(dmur)[2:5] = c('mu_1','mu_2','mu_3','mu_4')

ggplot(data.frame(dmur))+
  geom_line(aes(x,mu_1),col=5)+
  geom_line(aes(x,mu_2),col=3)+
  geom_line(aes(x,mu_3),col=2)+
  geom_line(aes(x,mu_4),col=4)+
  ggtitle('mu')
  
ggplot() + 
  geom_point(aes(x,Kr2))+
  ggtitle('Kr2')


dnr=cbind(x,nr)
colnames(dnr)[2:5] = c('nr_1','nr_2','nr_3','nr_4')

ggplot(data.frame(dnr))+
  geom_line(aes(x,nr_1),col=5)+
  geom_line(aes(x,nr_2),col=3)+
  geom_line(aes(x,nr_3),col='pink')+
  geom_line(aes(x,nr_4),col=4)+
  ggtitle('nr')



mean(tail(mur[!is.na(mur[,1]),1],500))
mean(tail(nr[!is.na(nr[,1]),1],500))
mean(tail(mur[!is.na(mur[,2]),2],500))
mean(tail(nr[!is.na(nr[,2]),2],500))
mean(tail(mur[!is.na(mur[,3]),3],500))
mean(tail(nr[!is.na(nr[,3]),3],500))

set.seed(416)
y_sample = c(rnorm(31, -1.07,0.3), rnorm(40,0.867,0.3), rnorm(14,-0.05,0.3))

plot(density(y_sample))
lines(density(Y))

plot(Kr2)

plot(nr[,1],type='l')
plot(mur[,1],type='l')
plot(nr[,2],type='l')
plot(mur[,2],type='l')
plot(nr[,3],type='l')
plot(mur[,3],type='l')
plot(nr[,4],type='l')
plot(mur[,4],type='l')

mean(mur[!is.na(mur[,3]),3])
mean(nr[!is.na(nr[,3]),3])

mean(mur[!is.na(mur[,4]),4])
mean(nr[!is.na(nr[,4]),4])

mean(mur[!is.na(mur[,5]),5])
mean(nr[!is.na(nr[,5]),5])

mean(mur[!is.na(mur[,6]),6])
mean(nr[!is.na(nr[,6]),6])

y_sample = c(rnorm(50,0.55,.1), rnorm(49,1.06,.1),
             rnorm(30,1.66,1),rnorm(86,-1.09,.1))
plot(density(y_sample))
plot(density(Y))
