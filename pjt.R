R = 3000 # Number of Gibbs iterations
tenit = R/10 # To be used later for the last 10% of the Gibbs iterations


# When using synthetic data, add the next 5 lines:
#-------------------------------------------------------------------------------------------
#mu0 = 2
#sigma0 = 5
#[Y, mugen] = synthetic(mu0, sigma0); % Function discussed and referenced in section 4.1.1
#Y=Y';
#Y=Y(1,:);
#%-------------------------------------------------------------------------------------------
#N = 300; % (Fixed) number of rows in all of the data sets
N = 200
Davis <- read.csv("Davis.csv")
#Y = Davis$weight
#Y = scale(Y)
#t = (expTest(:,1) - expTest(1,1)); % Equivalent to t = 0.1:0.1:30
t = seq(0.1,20,0.1)
#% Gaussian distribution
#fn = @(y,mu,sigma) 1/sqrt(2*pi*sigma)*exp(-0.5*(y-repmat(mu,1,size(y,2))).ˆ2/sigma);
#% Now ready to perform Gibbs sampling
#% Prior variance of mu
sigma0 = .1
#% Prior Dirichlet concentration of indicator probability p
c = 1e-3# % 1e-4 for Test 12
#% Likelihood variance of mu
sigma = .1
#% Initialization of the variables discussed
Kr = rep(1,R)
Kr2 = rep(1,R)
xr = matrix(rep(1,R*N),nrow=R)
mur = matrix(nrow=R,ncol=N)
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
    if(sum(i)==0){
      Kr2[r] = Kr[r]-1
      mu[k] = 10
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
Y = c(rnorm(50,1,.1), rnorm(50,0.5,.1),
      rnorm(50,-0.5,1),rnorm(50,-1,.1))
plot(Kr2)
plot(mur[,1],type='l')
plot(mur[,2],type='l')
plot(mur[,3],type='l')
plot(mur[,4],type='l')


