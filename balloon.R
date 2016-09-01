## Ballon example

# library(ftnonpar)

source('pen-s-functions.R')

data(balloon, package='ftnonpar')
y <- balloon
n <- length(y)
x <- (1:n)/n
plot(y~x, main='Balloon data', pch=19, col='gray', cex=.7)


# # ####### subsample of data ############
# x1 <- seq(1,4984,5)
# y1 <- y[x1] 
# n1 <- length(x1)
# x1 <- (1:n1)/n1 
# plot(y1~x1, main='Balloon data - subset 1', pch=19, col='gray', cex=.7)
# ##############

n1 <- 1500
set.seed(123)
x1 <- x[ii <- sort(sample(n, n1, repl=FALSE))]
y1 <- y[ii]
plot(y1~x1, main='Balloon data - subset 1', pch=19, col='gray', cex=.7)

x <- x1
y <- y1
n <- n1
p <- 3
num.knots <- max(5, min(floor(length(unique(x))/4),35))
knots <- quantile(unique(x),seq(0,1,length=num.knots+2))[-c(1,(num.knots+2))]
xpoly <- rep(1,n)
for (j in 1:p) xpoly <- cbind(xpoly,x^j)
xspline <- outer(x, knots, "-")
xspline <- pmax(xspline, 0)^p
X <- cbind(xpoly,xspline)
D <- diag(c(rep(0,ncol(xpoly)),rep(1,ncol(xspline))))

# N <- 500
# # space to return the selected lambdas
# opt.lam.m <- opt.lam.ls <- opt.lam.s <- rep(0, N)
# # space to return the vector of optimal GCV criteria
# gcvs.m <- gcvs.ls <- gcvs.s <- rep(0,N)
# # space to return the MSEs
# mse.m <- mse.ls <- mse.s <- rep(0,N)

# grid of lambdas to consider
lambdas <- seq(.001, 1, length=100)
ll <- length(lambdas)

# tuning constants for the S-estimator
cc <- 1.54764
b <- .5

# Cubic polynomials
NN <- 2 # ??

tmp.ls <- pen.ls.gcv(y, X, D, lambdas)
opt.lam.ls <- tmp.ls$lam
gcvs.ls <- tmp.ls$gcv

tmp.m <- pen.m.rcv(y=y, X=X, N=NN, D=D, lambdas=lambdas, num.knots=num.knots, p=p, epsilon=1e-6) 
opt.lam.m <- tmp.m$lam
gcvs.m <- tmp.m$gcv


tmp.s <- pen.s.gcv(y, X, D, lambdas, num.knots, p, NN, cc, b)
opt.lam.s <- tmp.s$lam
gcvs.s <- tmp.s$gcv

plot(y ~ x, main='Balloon data - subset 1', pch=19, col='gray', cex=.7, ylim=c(0,2.5))
lines(x, tmp.ls$yhat, lwd=2)
lines(x, tmp.s$yhat, col="blue", lwd=2)
lines(x, tmp.m$yhat, col="red", lwd=2)

plot(y ~ x, main='Balloon data - subset 1', pch=19, col='gray', cex=.7, ylim=c(0,2.5))
yhat.0 <- X %*% solve( t(X) %*% X + 1e-10 * D, t(X) %*% y)
lines(x, yhat.0, lwd=2, col='red')






### using spm function to estimate LS estimator####
library(SemiPar)
summary(m.plm<-spm(y~ f(x,lambda=0.1)))
#lines(m.plm,se=FALSE,col=3,lty=2)

