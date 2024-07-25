\name{interpolate}
\alias{interpolate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Resize array
}
\description{
Resize array
}
\usage{
interpolate(A, targetdim)
}
\arguments{
\item{A}{array or list}
\item{targetdim}{numeric; the size you want to resize}

}

\value{
the resized array
}

\author{
Gigi Wu
}

\examples{
im <- matrix(0, 558, 558)
im[186:372, ] <- 1
im[, 186:372] <- 1
#"im" is a matrix where the shape of nonzero element is cross.


##resize im
imRS = interpolate(im, c(32, 32))

dim(im)
dim(imRS)
}
