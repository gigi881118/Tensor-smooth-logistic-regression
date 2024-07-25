KTtoTensor <- function(x){
  if(class(x) != "list"){print("input must be list"); break}
  d = length(x)
  # print(paste0("d: ", d))
  rs = sapply(1:d, function(i) ncol(x[[i]]))
  if(all(rs == rs[1])){r = rs[1]}else{
    print("this is not a tensor."); break
  }
  # print(paste0("r: ", r))
  out = list()
  for (i in 1:r) {
    for (j in 1:d) {
      if(j == 1){out[[i]] = x[[j]][, i]}else{
        #out[[i]] = out[[i]] %o% x[[j]][, i]
        out[[i]] = outer(out[[i]], x[[j]][, i], FUN = "*")
      }
    }
  }
  beta_list = Reduce("+", out)
  return(beta_list)
}
