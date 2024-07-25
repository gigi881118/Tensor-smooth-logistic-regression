\name{KTtoTensor}
\alias{KTtoTensor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Turn the tensor in decomposition form to the regular tensor
}
\description{
Turn the tensor in decomposition form to the regular tensor (or multidimensional array)
}
\usage{
KTtoTensor(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the tensor in decompsition form which is a output of \code{TensorReg}}
}

\value{
a regular tensor in array form
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
}
