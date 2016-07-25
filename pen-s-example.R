
#############################################################################
#  Example script
#  R-code with functions implementing the penalised S-estimators used in the 
#  paper:  "S-Estimation for Penalized Regression Splines" 
#  published in the Journal of Computational and Graphical Statistics
#
#  Kukatharmini Tharmaratnam, Gerda Claeskens, Christophe Croux, and Matias 
#  Salibian-Barrera.
#  This file also includes a comparison with penalized least squares and 
#  penalized M-estimation methods.
#############################################################################

source('pen-s-functions.R')

set.seed(17)


n <- 100   # sample size
NN <- 500  # max. no. of iterations for S-estimator
NNN <- 500 # no. of initial candidates for the S-estimator
p0 <- .25   # prob of outlier generating mechanism
y0 <- 20   # location of outliers
sd0 <- 2   # sd of outliers

# grid of penalty parameters
# for the GCV search
lambdas <- seq(.1, 2, length=100)

# tuning constants for the S-estimator
cc <- 1.54764
b <- .5

# Cubic polynomials
p <- 3

# explanatory variables (x's)
x <- sort(runif(n,min=-1,max=1))

# true mean function
muf <- sin(pi*x)

# build the design matrix
num.knots <- max(5,min(floor(length(unique(x))/4),35))
knots <- quantile(unique(x),seq(0,1,length=num.knots+2))[-c(1,(num.knots+2))]
xpoly <-rep(1,n)
for (j in 1:p) xpoly<-cbind(xpoly,x^j)
xspline <- outer(x,knots,"-")
xspline <- pmax(xspline, 0)^p
#xspline <- xspline*(xspline>0)
X <- cbind(xpoly,xspline)

# penalty matrix
D <- diag(c(rep(0,ncol(xpoly)),rep(1,ncol(xspline))))

# error
e <- rnorm(n, mean=0, sd=.7)

# response variable
y <- muf + e

# outliers?
index <- (1:n)[ runif(n) < p0 ]
y[index] <- rnorm(length(index),y0,sd0)

# plot the data
plot(x,y, main='Comparison of penalized regression estimators')

# add the true mean function to the plot
lines(x, muf, lwd=3, col='black')

# penalized LS fit (with GCV)
tmp.ls <- pen.ls.gcv(y, X, D, lambdas)

# add the penalized LS estimate to the plot
lines(x, tmp.ls$yhat, lwd=3, col='red')

# penalized M-estimator (with robust CV)
tmp.m <- pen.m.rcv(y=y, X=X, N=NN, D=D, lambdas=lambdas, num.knots=num.knots, p=p, epsilon=1e-6) 

# add the penalized M estimate to the plot
lines(x, tmp.m$yhat, lwd=3, col='green')

# penalized S-estimator (with robust GCV)
tmp.s <- pen.s.rgcv(y=y, X=X, D=D, lambdas=lambdas, num.knots=num.knots, p=p, NN=NN, cc=cc, b=b, NNN=50)

# add the penalized S estimate to the plot
lines(x, tmp.s$yhat, lwd=3, col='blue')

legend(-.75, 15, legend=c("LS", "M", "S", "True"), 
col=c('red', 'green', 'blue', 'black'), lwd=3)


