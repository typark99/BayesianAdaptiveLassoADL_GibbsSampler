##################################################################################
## Run the Bayesian adaptive lasso for Autoregressive Distributed Lags model
## Taeyong Park, Feb 2016
##################################################################################



# MCMC LOOPING
for (M in 1:nIter)  {
  if (M %% 1000 == 0) print(M)
  
  #Full Conditional Posteriors  
  # beta
  dtau.inv   <- diag(1/tau.sq)
  cov.be     <- sig.sq * solve(XX+dtau.inv) 
  mean.be    <- 1/sig.sq * cov.be%*%t(X)%*%(Y.til)
  beta.p[M+1,] <- rmvnorm(1,mean=mean.be,sigma=cov.be)
  
  
  gam <- c()
  for (j in 1:nY){
    repeat{
      gam[j]  <- rinvGauss(1, nu=sqrt(lambda.post[M,j] * sig.sq/beta.p[M+1,j]^2), lambda=lambda.post[M,j])
      if (gam[j] > 0) break    	
    }
    tau.sqY[j] <- 1/gam[j]
    sh.lam      <- a + 1 
    sc.lam      <- 1/2*tau.sqY[j] + bY[j] 
    lambda.sqY[j] <- rgamma(1, shape=sh.lam, rate=sc.lam)
  }
  
  gam <- c()
  for (j in 1:nX){
    repeat{
      gam[j]  <- rinvGauss(1, nu=sqrt(lambda.post[M,j+nY] * sig.sq/beta.p[M+1,j+nY]^2), lambda=lambda.post[M,j+nY])
      if (gam[j] > 0) break    	
    }
    tau.sqX[j] <- 1/gam[j]
    sh.lam      <- a + 1 
    sc.lam      <- 1/2*tau.sqX[j] + bX[j] 
    lambda.sqX[j] <- rgamma(1, shape=sh.lam, rate=sc.lam)
  }
  
  gam <- c()
  if(nZ==0) {
    lambda.sqZ <- c() # When no non-dynamic variables are considered 
  } else {
    for (j in 1:nZ){
      repeat{
        gam[j]  <- rinvGauss(1, nu=sqrt(lambda.post[M,j+nY+nX] * sig.sq/beta.p[M+1,j+nY+nX]^2), lambda=lambda.post[M,j+nY+nX])
        if (gam[j] > 0) break    	
      }
      tau.sqZ[j] <- 1/gam[j]
      sh.lam      <- a + 1 
      sc.lam      <- 1/2*tau.sqZ[j] + bZ 
      lambda.sqZ[j] <- rgamma(1, shape=sh.lam, rate=sc.lam)
    }
  }

  tau.sq = c(tau.sqY, tau.sqX, tau.sqZ)
  lambda.sq = c(lambda.sqY, lambda.sqX, lambda.sqZ)
  tausq.post[M+1,] <- tau.sq 	
  lambda.post[M+1,] <- lambda.sq
  
  
  # sig.sq
  sh.sig     <- (n-1+nTheta)/2
  sc.sig     <- 1/2*t(Y.til-X%*%beta.p[M+1,])%*%(Y.til-X%*%beta.p[M+1,])+ 1/2*t(beta.p[M+1,])%*%diag(1/lambda.sq)%*%beta.p[M+1,]
  sig.sq     <- rinvgamma(1, shape=sh.sig, scale=sc.sig)
  sigsq.post <- c(sigsq.post, sig.sq)
  
}