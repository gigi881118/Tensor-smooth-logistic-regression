\name{TensorRegNet}
\alias{TensorRegNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tensor Regression with Lasso or Smooth Regularization
}
\description{
Fit a generalized tensor model with Lasso or smooth regularization.
\code{TensorRegNet} is used to fit linear and logistic tensor regression models with regularization.
}
\usage{
TensorRegNet = function(
  X,
  M,
  y,
  r,
  dist,
  lambda,
  smooth = FALSE,
  penparam = 1,
  B0 = NULL,
  Display = FALSE,
  BurninMaxIter = 20,
  BurninTolFun = 1e-2,
  BurninReplicates = 5,
  PenaltyMaxIter = 100,
  PenaltyTolFun = 1e-3,
  warn = FALSE,
  wts  = NULL
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

\item{lambda}{
penalty tuning constant}

\item{smooth}{
logical; smooth tensor regression coefficient or not}

\item{penparam}{
the index parameter for elastic net mixing parameter \eqn{\alpha}, \eqn{\alpha in [0, 1]}}

\item{B0}{
starting point, it must be the tensor in decompsition form which is a output of \code{TensorReg}
}

\item{Display}{
logical; whether to show the details of iteration}

\item{BurninMaxIter}{
maximum iteration for the burn-in runs, default is 20}

\item{BurninTolFun}{
tolerance for the burn-in runs, default is 1e-2}

\item{BurninReplicates}{
number of the burn-in runs, default is 5}

\item{PenaltyMaxIter}{
maximum iteration at penalization stage, default is 100}

\item{PenaltyTolFun}{
tolerenceat penalization stage, default is 1e-3}

\item{warn}{
logical; whether to show the warning during the iterations}

\item{wts}{
observation weights; default is ones for each obs}
}

\value{

\item{beta0_final}{
regression coefficients for the auxiliary covariates
}
\item{beta_final}{
a tensor of regression coefficients for array variates in decompsition form
}
\item{beta_scale}{
a tensor of the scaling constants for the array coefficients}

\item{glmstats}{
GLM regression summary statistics from last iter.
}

}
\references{
Feng, L., Bi, X., & Zhang, H. (2021). Brain regions identified as being associated with verbal reasoning through the use of imaging regression via internal variation. Journal of the American Statistical Association, 116(533), 144-158.

Zhou, H., Li, L., & Zhu, H. (2013). Tensor regression with applications in neuroimaging data analysis. Journal of the American Statistical Association, 108(502), 540-552.
}
\author{
Gigi Wu
}

\examples{
load('./example/im.rdata')
#"im" is a matrix where the shape of nonzero element is cross.

## IM: Resize
imRS = interpolate(im, c(32, 32))
b = matrix(0, 64, 64)
b[(dim(b)[1]/4):((dim(b)[1]/4)+dim(imRS)[1]-1),
  (dim(b)[2]/4):((dim(b)[2]/4)+dim(imRS)[2]-1)] <- imRS


#---------- ex1: y is cont. -----------
## parameter setting
p0 = 5
b0 = matrix(rep(1, p0), p0, 1)
p1 = dim(b)[1] ; p2 = dim(b)[2]
n = 500 #sample size
r = 2
distribution = gaussian()

## variates setting
X = array(rnorm(n*p0), dim=c(n, p0)) #n-by-p0 regular design matrix
M = array(rnorm(p1*p2*n), dim=c(p1,p2,n)) #p1-by-p2-by-p3-by-n tensor variates
mu = X \%*\% b0 + t(t(matrix(b)) \%*\% matrix(M, ncol = n))
sigma = 1
y = mu + sigma*array(rnorm(n), dim=c(n, 1))

## Run
#TR
TR_cont = TensorReg(X,M,y,r, dist = distribution)
beta_TR_cont = TR_cont$beta_final
B0reg = beta_TR_cont #initial B0

#TR LASSO
##grid search lambda
lambdas = c(0.01, 0.03, 0.05, 0.1, 0.5, 1)
Lambda_r1 = sapply(lambdas, function(L){
  TRsp <- TensorRegNet(X, M, y, r, dist = distribution,
                             lambda = L, B0 = B0reg, smooth = F)
  return(TRsp$glmstats$BIC)
})
(finalLambda_r1 = lambdas[which.min(Lambda_r1)])

TRLASSO_cont = TensorRegNet(X,M,y,r, dist = distribution,
                                  lambda = finalLambda_r1, B0 = B0reg,
                                  smooth = F)
beta_TRLASSO_cont = TRLASSO_cont$beta_final


#TR Smooth
lambdas2 = c(0.05, 0.5, 1, 5, 10, 30)
Lambda_r2 = sapply(lambdas2, function(L){
  TRsp <- TensorRegNet(X, M, y, r, dist = distribution,
                             lambda = L, B0 = B0reg, smooth = T)
  return(TRsp$glmstats$BIC)
})
(finalLambda_r2 = lambdas2[which.min(Lambda_r2)])

TRsmooth_cont = TensorRegNet(X,M,y,r, dist = distribution,
                                   lambda = finalLambda_r2, B0 = B0reg,
                                   smooth = T)
beta_TRsmooth_cont = TRsmooth_cont$beta_final


## plot the coef.
windows()
par(mfrow=c(2, 2))

plot(as.cimg(abs(b - 1)), xlim = c(0, 64), ylim = c(0, 64),
     main = "True Signal")

plot(as.cimg(abs(KTtoTensor(beta_TR_cont) - 1)), xlim = c(0, 64), ylim = c(0, 64),
     main = "TR"
)

plot(as.cimg(abs(KTtoTensor(beta_TRLASSO_cont) - 1)), xlim = c(0, 64), ylim = c(0, 64),
     main = "TR LASSO"
)

plot(as.cimg(abs(KTtoTensor(beta_TRsmooth_cont) - 1)), xlim = c(0, 64), ylim = c(0, 64),
     main = "TR Smooth"
)


#---------- ex2: y is bin. -----------
## parameter setting
p0 = 5
b0 = matrix(rep(1, p0), p0, 1)
n = 500 #sample size
r = 2
distribution = binomial()

## variates setting
X = array(rnorm(n*p0), dim=c(n, p0)) #n-by-p0 regular design matrix
M = array(rnorm(p1*p2*n), dim=c(p1,p2,n)) #p1-by-p2-by-n matrix variates
mu = X \%*\% b0 + t(t(matrix(b)) \%*\% matrix(M, ncol = n))
sigma = 1
y = as.matrix(sapply(mu, function(x) rbinom(1, 1, 1/(1+exp(-x)))))

## Run
#TR
TR_bin = TensorReg(X,M,y,r, dist = distribution)
beta_TR_bin = TR_bin$beta_final
B0reg = beta_TR_bin #initial B0

#TR LASSO
##grid search lambda
lambdas = c(0.01, 0.03, 0.05, 0.5, 1)
Lambda_r1 = sapply(lambdas, function(L){
  TRsp <- TensorRegNet(X, M, y, r, dist = distribution,
                             lambda = L, B0 = B0reg)
  return(TRsp$glmstats$BIC)
})
(finalLambda_r1 = lambdas[which.min(Lambda_r1)])
TRLASSO_bin = TensorRegNet(X,M,y,r, dist = distribution,
                                 lambda = finalLambda_r1, B0 = B0reg,
                                 smooth = F)
beta_TRLASSO_bin = TRLASSO_bin$beta_final



#TR Smooth
##grid search lambda
lambdas2 = c(0.05, 0.5, 1, 5, 10, 30)
Lambda_r2 = sapply(lambdas2, function(L){
  TRsp <- TensorRegNet(X, M, y, r, dist = distribution,
                             lambda = L, B0 = B0reg, smooth = T)
  return(TRsp$glmstats$BIC)
})
(finalLambda_r2 = lambdas[which.min(Lambda_r2)])
TRsmooth_bin = TensorRegNet(X,M,y,r, dist = distribution,
                                  lambda = finalLambda_r2, B0 = B0reg,
                                  smooth = T)
beta_TRsmooth_bin = TRsmooth_bin$beta_final


#plot the coef.
windows()
par(mfrow=c(2, 2))

plot(as.cimg(abs(b - 1)), xlim = c(0, 64), ylim = c(0, 64),
     main = "True Signal")

plot(as.cimg(abs(KTtoTensor(beta_TR_bin) - 1)), xlim = c(0, 64), ylim = c(0, 64),
     main = "TR"
)

plot(as.cimg(abs(KTtoTensor(beta_TRLASSO_bin) - 1)), xlim = c(0, 64), ylim = c(0, 64),
     main = "TR LASSO"
)

plot(as.cimg(abs(KTtoTensor(beta_TRsmooth_bin) - 1)), xlim = c(0, 64), ylim = c(0, 64),
     main = "TR Smooth"
)



}

