interpolate <- function(A, targetdim){

  # turn array into tensor structure
  p = c()
  if(any(class(A) == "array")){
    p = dim(A)
  }else{
    if(class(A) == "list"){
      for (i in 1:length(A)) {
        p = c(p, nrow(A[[i]]))
      }
    }
  }
  D = length(p)

  # argument checking
  if(length(targetdim) != D) stop('targetdim does not match dimension of A')

  U = list()
  for (d in 1:D) {
    xi = seq(1, p[d], length = targetdim[d])
    j1 = floor(xi)
    j2 = j1 + 1
    w1 = j2 - xi
    w2 = 1 - w1
    j2[length(j2)] = p[d]
    #the interpolation matrix: targetdim(d)-by-p(d)
    U[[d]] = sparseMatrix(i = rep(1:targetdim[d], times=2),
                          j = c(j1, j2),
                          x = c(w1, w2),
                          dims = c(targetdim[d], p[d]))
  }


  if(any(class(A) == "array")){
    B = new("Tensor",D, p, data=A)
    for (d in 1:D) {
      B = ttm(B, as.matrix(U[[d]]), m=d)
    }
    #if(is.integer(A)){
    #  return(round(B@data))
    #}else{
    #  return(B@data)
    #}
    return(B@data)
  }else{
    if(is.list(A)){
      B = A
      for (d in 1:D) {
        #B[[d]] = as.matrix(U[[d]])%*%B[[d]]# resize factor matrix
        B[[d]] = tcrossprod(as.matrix(U[[d]]), t(B[[d]]))
      }
      #if(is.integer(A)){
      #  return(round(B))
      #}else{
      #  return(B)
      #}
      return(B)
    }
  }
}



