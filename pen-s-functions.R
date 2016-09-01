
################################################################################
#  R-code with functions implementing the penalised S-estimators used in the 
#  paper:  "S-Estimation for Penalized Regression Splines" 
#  published in the Journal of Computational and Graphical Statistics
#
#  Kukatharmini Tharmaratnam, Gerda Claeskens, Christophe Croux, and Matias 
#  Salibian-Barrera.
#  This file also includes a comparison with penalized least squares and 
#  penalized M-estimation methods.
#  An example script is included at the end of this file.
################################################################################
#
#
# Define the function for Penalized S-estimators
pen.s <- function(y, X, N, D, lambda, num.knots, p, beta1,
                  Sbeta1, cc=1.54764, b=.5, epsilon=1e-6) 
{ 
  betahats <- matrix(ncol=num.knots+2+p-1+1,nrow=N) 
  #store the values in a matrix 
  betahats[1,] <- c(beta1,Sbeta1)
  beta <- beta1 
  for (i in 2:N)
  { 
    #update Sbeta conditional on beta 
    r <- as.vector(y-X%*%beta)
    Sbeta <- s.scale(r, cc=cc, b=b, N, ep=1e-4)
    rs <- r / Sbeta
    Wbeta <- Psi(rs, cc) / rs
    taubeta <- n*(Sbeta)^2 / sum( r^2 * Wbeta ) 
    #update beta conditional on Sbeta from above
    beta <- solve(t(X*Wbeta)%*%X +(D*lambda/taubeta))%*% t(X*Wbeta)%*% y
    betahats[i,] <- c(beta,Sbeta)
    ifelse(((norm(betahats[i,]-betahats[i-1,])/norm(betahats[i-1,]))<epsilon),break,next)
  }
  return(list(outmatrix=betahats[1:i,1:(num.knots+2+p-1+1)],
              estimates=betahats[i,1:(num.knots+2+p-1)],
              scale=betahats[i,num.knots+2+p-1+1],
              iterations=i,weights=Wbeta))
}

# Robust GCV search for penalized S-estimators 
pen.s.rgcv <- function(y, X, D, lambdas, num.knots, p, NN, cc, b,NNN) 
{ 
  ll <- length(lambdas) 
  # GCVs for the S estimator 
  best.gcv <- +Inf 
  rgcvs <- rep(0,ll) 
  for(i in 1:ll) 
  {
    uu<- initial.S(y, X, D, lambdas[i], num.knots, p, NN, cc, b, NNN) 
    rgcvs[i] <- rgcv(uu, y, X, D, lambdas[i]) 
    if( rgcvs[i] <= best.gcv ) 
    {
      best.uu <- uu
      best.gcv <- rgcvs[i]
    }
  }
  # find the best lambda
  #rlam <- best.gcv 
  rlam<- lambdas[rgcvs==best.gcv]
  yhat.s <- as.vector( X %*% best.uu$estimates ) 
  return(list(yhat = yhat.s, lam=rlam, gcv=min(rgcvs), 
              iter.s=best.uu$iterations))
}

# Define the penalized S-estimators with different initial candidates
initial.S<- function(y, X, D,lambda, num.knots, p, NN, cc, b, NNN) 
{
  uubeta <- matrix(0, ncol=(NNN+2),nrow=num.knots+2+p-1)  
  # To get best beta w.r.t objective function 
  uuiteration <- rep(0,(NNN+2)) # Needed to use pen.s.rgcv instead of pen.s() 
  uuscale <- rep(0,(NNN+2))     
  uuweights <- matrix(0, ncol=(NNN+2),nrow=n) 
  objval <- rep(0,(NNN+2))              # To get min of objval 
  # Initial candidates from Resampling 
  for (ii in 1:NNN)
  {
    indices <- sample(n,num.knots+2+p-1+1)
    Xs <- X[indices,]
    ys <- y[indices] 
    init <- pen.ls(ys, Xs, D, lambda)
    uu1 <- pen.s(y,X,20,D,lambda,num.knots,p,init$beta,init$Sbeta,cc=cc,b=b)
    uubeta[,ii]<- as.vector(uu1$estimates)
    uuscale[ii] <- uu1$scale
    uuweights[,ii] <- as.vector(uu1$weights)
    uuiteration[ii] <- uu1$iterations
    objval[ii] <- ((n*(uu1$scale^2))+(lambda*as.numeric(t(uu1$estimates)%*% 
                                                          D %*% uu1$estimates)))
  }
  # Initial candidates from M-estimator
  initM <- pen.m(y, X, N=NN, D, lambda, num.knots, p) 
  uuM <-pen.s(y,X,20,D,lambda,num.knots,p,initM$outmbeta,initM$sigma,cc=cc,b=b)
  uubeta[,(NNN+1)]<- as.vector(uuM$estimates)
  uuscale[(NNN+1)] <- uuM$scale
  uuweights[,(NNN+1)] <- as.vector(uuM$weights)
  uuiteration[(NNN+1)] <- uuM$iterations
  objval[(NNN+1)] <- ((n*(uuM$scale^2))+(lambda* as.numeric(t(uuM$estimates) 
                                                            %*% D %*% uuM$estimates)))
  # Initial candidates from LS-estimator
  initLS <- pen.ls(y, X, D, lambda)
  uuLS <-pen.s(y,X,20,D,lambda,num.knots,p,initLS$beta,initLS$Sbeta,cc=cc,b=b)
  uubeta[,(NNN+2)]<- as.vector(uuLS$estimates)
  uuscale[(NNN+2)] <- uuLS$scale
  uuweights[,(NNN+2)] <- as.vector(uuLS$weights)
  uuiteration[(NNN+2)] <- uuLS$iterations
  objval[(NNN+2)] <-((n*(uuLS$scale^2))+(lambda*as.numeric(t(uuLS$estimates) 
                                                           %*% D %*% uuLS$estimates)))
  # find the best estimators with respect to objective function
  bestbeta <-as.vector( uubeta[ ,objval == min(objval) ])
  weights <-as.vector(uuweights[ ,objval == min(objval) ])
  scale <- uuscale[objval == min(objval)]
  iterations <- uuiteration[objval == min(objval)]
  return(list(estimates=bestbeta, scale=scale, weights=weights, 
              iterations=iterations))
}

#Define Scale function  
s.scale <- function(r, cc=1.54764, b=.5, max.it=1000, ep=1e-4) 
{
  s1 <- mad(r)
  if(abs(s1)<1e-10) return(s1)
  s0 <- s1 + 1
  it <- 0
  while( ( abs(s0-s1) > ep ) && (it < max.it) ) {
    it <- it + 1
    s0 <- s1
    s1 <- s0*mean(Rho(r/s0,cc=cc))/b
  }
  return(s1)
}

# Define robust generalized cross validation function for S-penalized regression
rgcv <- function(uu, y, X, D, lambda) 
{ 
  # uu has the fit returned by pen.s()
  # y is the response vector
  # X is the big design matrix
  # D is the penalty matrix
  # lambda is the value of the penalty constant to be evaluated
  n <- length(y)
  nw <- sum( uu$weights > 0 )
  r <- as.vector(y - X %*% uu$estimates )
  aa <- n * uu$scale^2 / sum( r^2 * uu$weights )
  sw <- sqrt(uu$weights)
  H <- (X*sw)%*%solve(t(X*uu$weights)%*%X+lambda/aa*D)%*%t(X*sw)
  return(nw * sum( r^2 * uu$weights ) / (nw - sum(diag(H)))^2 ) 
}


#Define rho function
Rho<- function(x, cc)
{
  U <- x/cc
  U1 <- 3 * U^2 - 3 * U^4 + U^6
  U1[abs(U) > 1] <- 1
  return(U1)
}

#Define psi function
Psi<-function(x, cc) 
{
  U <- x/cc
  U1 <- 6/cc * U * (1 - U^2)^2
  U1[abs(U) > 1] <- 0
  return(U1)
}

#Define norm function  
norm <- function(a) sqrt(sum(a^2))


######
# Comparison with penalized least squares estimation.
######

# Define the function for Penalized LS-estimators
pen.ls <- function(y, X, D, lambda)
{
  beta.ls <- as.vector(solve( t(X) %*% X + lambda * D ) %*% t(X) %*% y )
  Sbeta.ls <- mad( y - X %*% beta.ls)
  return(list(beta=beta.ls,Sbeta=Sbeta.ls))
}

# Define generalized cross validation function for LS-penalized regression
gcv <- function(y, X, D, lambda) 
{ 
  # y is the response vector
  # X is the big design matrix
  # D is the penalty matrix
  # lambda is the value of the penalty constant to be evaluated
  tmp <- solve( t(X) %*% X + lambda * D, t(X)) 
  beta <- as.vector( tmp %*% y )
  n <- length(y)
  r <- as.vector(y - X %*% beta)
  H <- X %*% tmp
  return( n * sum( r^2 ) / (n - sum(diag(H)))^2 )
}

# GCV search for penalized LS-estimators
pen.ls.gcv <- function(y, X, D, lambdas) 
{
  ll <- length(lambdas)
  # GCVs for the LS estimator
  gcvs <- rep(0, ll)
  for(i in 1:ll) 
  {
    gcvs[i] <- gcv(y, X, D, lambdas[i])
  }
  # find the best lambda
  lam <- max( lambdas[ gcvs == min(gcvs) ] )
  beta.ls <- as.vector(solve( t(X) %*% X + lam * D, t(X) %*% y ))
  # get the LS estimated mean
  yhat.ls <- as.vector( X %*% beta.ls )
  Sbeta.ls <- mad( y - yhat.ls)
  return(list(beta=beta.ls, Sbeta=Sbeta.ls, yhat = yhat.ls, lam=lam, gcv=min(gcvs)))
}

######
# Comparison with penalized M-estimation.
######

# Define the function for penalized M-estimators for fixed lambda - (Proposed by # Oh-Lee 2007)
pen.m<- function(y, X, N, D, lambda, num.knots, p, epsilon=1e-6) 
{
  results <- matrix(ncol=n+1,nrow=N) #store the values in matrix
  # start with penalized LS
  tmp <- pen.ls(y, X, D, lambda)
  beta1 <- as.vector( tmp$beta )
  mhat1 <- as.vector( X %*% tmp$beta )
  sigma1 <- tmp$Sbeta
  results[1,] <- c(mhat1, sigma1)
  mhat <- mhat1
  mbetaresults <- matrix(ncol=num.knots+2+p-1, nrow=N)
  mbetaresults[1,] <- c(beta1)
  mbeta <- beta1
  for (j in 2:N)
  { 
    res <- as.vector(y-X%*%mbeta)
    # sigma <- 1.4826*median(abs(res))
    sigma <- mad(res)
    cval <- 1.345*sigma
    psi1 <- ifelse( abs(res)<=cval, 2*res, 2*cval*sign(res) )
    z <- mhat + (psi1/2)
    mbeta <- solve( t(X)%*%X + D*lambda ) %*% t(X) %*% z
    mhat <- as.vector( X %*% solve( t(X)%*%X + D*lambda ) %*% t(X) %*% z )
    results[j,] <- c(mhat,sigma)
    mbetaresults[j,] <- c(mbeta)
    ifelse(((norm(mbetaresults[j,]-mbetaresults[j-1,])/norm(mbetaresults[j-
                                                                           1,]))<epsilon),break,next)
  }
  return(list(outmbeta=as.vector(mbetaresults[j,]),
              sigma=as.vector(results[j,n+1]), iterations=j))
}


# Define robust cross validation function for M-penalized regression-
# Cantoni and Ronchetti (2001)
mrcv <- function(mm, y, X, D, lambda,n) 
{ 
  # mm has the fit returned by pen.m()
  # y is the response vector
  # X is the big design matrix
  # D is the penalty matrix
  # lambda is the value of the penalty constant to be evaluated
  #n <- length(y)
  res <- as.vector(y - X %*% mm$outmbeta )
  sigma <- mad(res)
  cval <- 1.345*sigma
  psi1 <- ifelse( abs(res)<=cval, 2*res, 2*cval*sign(res) )
  psi1dash <- ifelse( abs(res)<=cval, 2,0 )
  Epsi1dash <- sum(psi1dash)/n
  II<- diag(c(rep(1,ncol(X)))) 
  SS <-X %*% solve(II+ lambda * (sigma/Epsi1dash)* D, t(X))
  return( mrcv=1/n * (sigma^2/Epsi1dash^2)* sum(psi1^2/ (1- diag(SS))^2 )) 
}

# Robust CV search for penalized M-estimators 
pen.m.rcv <- function(y, X, NN, D, lambdas, num.knots, p, epsilon=1e-6) 
{ 
  ll <- length(lambdas) 
  # MCVs for the M estimator 
  best.cv <- +Inf 
  mrcvs <- rep(0,ll) 
  for(i in 1:ll) 
  {
    mm<- pen.m(y, X, NN, D, lambdas[i], num.knots, p, epsilon=1e-6) 
    mrcvs[i] <- mrcv(mm, y, X, D, lambdas[i],n) 
    if( mrcvs[i] <= best.cv ) 
    {
      best.mm <- mm
      best.cv <- mrcvs[i]
    }
  }
  # find the best lambda
  #rlam <- best.gcv 
  rlam <- max( lambdas[mrcvs==best.cv] )
  yhat.m <- as.vector( X %*% best.mm$outmbeta ) 
  return(list(yhat = yhat.m, lam=rlam, gcv=min(mrcvs),
              sigma.m=best.mm$sigma, iterations=best.mm$iterations))
}

