S-estimation for penalized regression splines
================
Matias Salibian
2016-08-31

S-estimators for penalized regression splines
---------------------------------------------

This repository contains code implementing the penalised S-regression estimators as proposed in "S-estimation for penalised regression splines" [Tharmaratnam et al (2010)](http://dx.doi.org/10.1198/jcgs.2010.08149).

All the necessary functions are in the file `pen-s-functions.R`. The following example illustrates how the code is used.

``` r
source("pen-s-functions.R")

set.seed(17)


n <- 100  # sample size
NN <- 500  # max. no. of iterations for S-estimator
NNN <- 500  # no. of initial candidates for the S-estimator
p0 <- 0.25  # prob of outlier generating mechanism
y0 <- 20  # location of outliers
sd0 <- 2  # sd of outliers

# grid of penalty parameters for the GCV search
lambdas <- seq(0.1, 2, length = 100)

# tuning constants for the S-estimator
cc <- 1.54764
b <- 0.5

# Cubic polynomials
p <- 3

# explanatory variables (x's)
x <- sort(runif(n, min = -1, max = 1))

# true mean function
muf <- sin(pi * x)

# build the design matrix
num.knots <- max(5, min(floor(length(unique(x))/4), 35))
knots <- quantile(unique(x), seq(0, 1, length = num.knots + 2))[-c(1, (num.knots + 
    2))]
xpoly <- rep(1, n)
for (j in 1:p) xpoly <- cbind(xpoly, x^j)
xspline <- outer(x, knots, "-")
xspline <- pmax(xspline, 0)^p
# xspline <- xspline*(xspline>0)
X <- cbind(xpoly, xspline)

# penalty matrix
D <- diag(c(rep(0, ncol(xpoly)), rep(1, ncol(xspline))))

# error
e <- rnorm(n, mean = 0, sd = 0.7)

# response variable
y <- muf + e

# outliers?
index <- (1:n)[runif(n) < p0]
y[index] <- rnorm(length(index), y0, sd0)

# plot the data
plot(x, y, main = "Comparison of penalized regression estimators")

# add the true mean function to the plot
lines(x, muf, lwd = 3, col = "black")

# penalized LS fit (with GCV)
tmp.ls <- pen.ls.gcv(y, X, D, lambdas)

# add the penalized LS estimate to the plot
lines(x, tmp.ls$yhat, lwd = 3, col = "red")

# penalized M-estimator (with robust CV)
tmp.m <- pen.m.rcv(y = y, X = X, N = NN, D = D, lambdas = lambdas, num.knots = num.knots, 
    p = p, epsilon = 1e-06)

# add the penalized M estimate to the plot
lines(x, tmp.m$yhat, lwd = 3, col = "green")

# penalized S-estimator (with robust GCV)
tmp.s <- pen.s.rgcv(y = y, X = X, D = D, lambdas = lambdas, num.knots = num.knots, 
    p = p, NN = NN, cc = cc, b = b, NNN = 50)

# add the penalized S estimate to the plot
lines(x, tmp.s$yhat, lwd = 3, col = "blue")

legend(-0.75, 15, legend = c("LS", "M", "S", "True"), col = c("red", "green", 
    "blue", "black"), lwd = 3)
```

![](README_files/figure-markdown_github/example-1.png)
