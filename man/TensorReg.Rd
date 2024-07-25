\name{TensorReg}
\alias{TensorReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tensor Regression
}
\description{
Fit generalized tensor regression models using the auxiliary covariate matrix X,  multidimensional array, that is tensor, covariates M and cont/binary response y.
\code{TensorReg} is used to fit linear and logistic tensor regression models.
}
\usage{
TensorReg(
  X = NULL,
  M,
  y,
  r,
  dist,
  B0 = NULL,
  Display = FALSE,
  MaxIter = 100,
  Replicates = 5,
  TolFun = 1e-4,
  wts = rep(1, n)
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{
NULL; n-by-p0 auxiliary covariate matrix}

\item{M}{
array variates (or tensors) with dim(M) = [p1,p2,...,pd,n]}

\item{y}{
n-by-1 respsonse vector;
if the respsonse is binary, y is coded in 0/1}

\item{r}{
rank; positive integer}

\item{dist}{
distribution; '\code{gaussian()}' or '\code{binomial()}'}

\item{B0}{
starting point, it must be the tensor in decompsition form which is a output of \code{TensorReg}}

\item{Display}{
logical; whether to show the details of iteration}

\item{MaxIter}{
maximum iteration, default is 100
}

\item{Replicates}{
number of intitial points to try, default is 5
}

\item{TolFun}{
tolerence in objective value, default is 1e-4
}

\item{wts}{
observation weights; default is ones for each obs
}
}

\value{
\item{beta0_final}{
regression coefficients for the auxiliary covariates
}
\item{beta_final}{
a tensor of regression coefficients for array variates in decompsition form
}
\item{glmstats_final}{
GLM regression summary statistics from last iter.
}

\item{dev_final}{
deviance of final model
}

%% ...
}
\references{
Zhou, H., Li, L., & Zhu, H. (2013). Tensor regression with applications in neuroimaging data analysis. Journal of the American Statistical Association, 108(502), 540-552.
}
\author{
Gigi Wu
}
\examples{
load('./example/im.rdata') #"im" is a matrix where the shape of nonzero element is cross.

##resize im
imRS = interpolate(im, c(32, 32))
b = matrix(0, 64, 64)
b[(dim(b)[1]/4):((dim(b)[1]/4)+dim(imRS)[1]-1),
  (dim(b)[2]/4):((dim(b)[2]/4)+dim(imRS)[2]-1)] <- imRS

#--------- ex1: y is cont. -------
## parameter setting
p0 = 5
p1 = dim(b)[1] ; p2 = dim(b)[2]
b0 = matrix(rep(1, p0), p0, 1)
n = 500 #sample size
r = 2
distribution = gaussian()

## variates setting
X = array(rnorm(n*p0), dim=c(n, p0)) #n-by-p0 regular design matrix
M = array(rnorm(p1*p2*n), dim=c(p1,p2,n)) #p1-by-p2-by-p3-by-n tensor variates
mu = X \%*\% b0 + t(t(matrix(b)) \%*\% matrix(M, ncol = n))
sigma = 1
y = mu + sigma*array(rnorm(n), dim=c(n, 1))

## run
TR_cont = TensorReg(X,M,y,r, dist = distribution)
beta_TR_cont = TR_cont$beta_final #tensor covariates in decompositon form
beta = KTtoTensor(beta_TR_cont) #tensor covariates in array form


## plot the coef.
windows()
par(mfrow=c(1, 2))

plot(as.cimg(abs(b - 1)), xlim = c(0, 64), ylim = c(0, 64),
     main = "True Signal")

plot(as.cimg(abs(KTtoTensor(beta_TR_cont) - 1)), xlim = c(0, 64), ylim = c(0, 64),
     main = "TR"
)

#--------- ex2: y is bin. -------
## parameter setting
p0 = 5
p1 = dim(b)[1] ; p2 = dim(b)[2]
b0 = matrix(rep(1, p0), p0, 1)
n = 500 #sample size
r = 2distribution = binomial()

## variates setting
X = array(rnorm(n*p0), dim=c(n, p0)) #n-by-p0 regular design matrix
M = array(rnorm(p1*p2*n), dim=c(p1,p2,n)) #p1-by-p2-by-n matrix variates
mu = X \%*\% b0 + t(t(matrix(b)) \%*\% matrix(M, ncol = n))
sigma = 1
y = as.matrix(sapply(mu, function(x) rbinom(1, 1, 1/(1+exp(-x)))))

## Run
TR_bin = TensorReg(X,M,y,r, dist = distribution)
beta_TR_bin = TR_bin$beta_final

## plot the coef.
windows()
par(mfrow=c(1, 2))

plot(as.cimg(abs(b - 1)), xlim = c(0, 64), ylim = c(0, 64),
     main = "True Signal")

plot(as.cimg(abs(KTtoTensor(beta_TR_bin) - 1)), xlim = c(0, 64), ylim = c(0, 64),
     main = "TR"
)
