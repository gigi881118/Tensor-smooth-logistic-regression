TensorReg = function(X=NULL,M,y,r,dist,
                       B0 = NULL,
                       Display = FALSE, #'iter'
                       MaxIter = 100,
                       Replicates = 5,
                       TolFun = 1e-4,
                       wts = rep(1, n)){
  #r=1;dist=gaussian(); B0 = NULL; Display = 'off'; #'iter'
  #MaxIter = 100;Replicates = 5;  TolFun = 1e-4; wts = rep(1, n)

  #----------- check dimensions ---------------
  if(is.null(X)){
    dimX = c(tail(dim(M), 1), 1)
    X = array(rep(1, prod(dimX)), dim = dimX)
  }

  n <- dim(X)[1]
  p0 <- dim(X)[2]
  d = length(dim(M)) - 1      # dimension of array variates
  p = dim(M)

  if (p[length(p)] != n) {
    stop("dimension of M does not match that of X!", call. = FALSE)
  }

  if (n < p0 || n < r * max(p[-length(p)])) {
    stop("sample size n is not large enough to estimate all parameters!", call. = FALSE)
  }

  #------------ mode-d matricization of TM ---------------
  TM = M

  mem_info <- gc()
  total_mem <- as.numeric(memory.limit()) * 1024 * 1024 #turn megabytes to bytes
  available_mem <- total_mem - mem_info[2]

  # CAUTION: may cause out of memory on Linux
  if(d*(8*prod(dim(TM)))<.75*available_mem){
    Md <- list()
    for (dd in 1:d) {
      Md[[dd]] = tenmat(TM, c(d+1, dd))
    }
  }

  #------- check user-supplied initial point: B0 ----------
  if(!is.null(B0)){
    # ignore requested multiple initial points
    #Replicates = 1
    #print("Replicates = 1")

    # check dimension
    px = c()
    if(any(class(B0) == "array")){
      px = dim(B0)
    }else{
      if(class(B0) == "list"){
        for (i in 1:length(B0)) {
          px = c(px, nrow(B0[[i]]))
        }
      }
    }
    ndimB = length(px)
    if(ndimB != d) stop('dimension of B0 does not match that of data!')

    # turn B0 into a tensor (if it's not)
    if(class(B0) == "numeric") B0 = as.tensor(B0)

    #resize to compatible dimension (if it's not)
    if(sum(px != p[-length(p)]) > 0){
      B0 = interpolate(B0, p[-length(p)]) ##edit
    }

    # perform CP decomposition if it's not a ktensor of correct rank
    if(any(class(B0) == "array") |
       (class(B0) == "list" & ncol(B0[[1]]) != r)) {
      B0 = cp(as.tensor(B0), num_components = r)
    }
  }

  #--- turn off warnings from glmfit_priv ------------

  #Available parameter name/value pairs are:
  #'B0': starting point, it can be a numeric array or a tensor
  #'Display': 'off' <default> or 'iter'
  #'MaxIter': maximum iteration, default is 100
  #'Replicates': # of intitial points to try, default is 5
  #'TolFun': tolerence in objective value, default is 1e-4
  #'weights': observation weights, default is ones for each obs.
  #.....................................................
  #B0 = NULL
  #Display = 'off' #'iter'
  #MaxIter = 100
  #Replicates = 5
  #TolFun = 1e-4
  #wts = rep(1, n)
  #................................................

  #--------- pre-allocate variables------------
  glmstats = list()
  dev_final = Inf


  #------ loop for various intial points -----------
  for(rep in 1:Replicates){
    if(!is.null(B0)){
      beta = B0
    }else{
      # initialize tensor regression coefficients from uniform (-1,1)
      #.............................................................
      ## Initialize: (alpha_0, gamma_0 = argmax...) <Algorithm_1>
      #beta <- lapply(1:d, function(j) matrix(1 - 2*runif(p[j]*r), ncol = r))
      beta <- lapply(1:d, function(j) matrix(runif(p[j]*r, -1, 1), ncol = r))
      beta.lambda <- rep(1, r) ##weight for each rank is same
    }

    # main loop
    for (iter in 1:MaxIter) {
      ## Repeat: <Algorithm_1>
      #.............................................................

      # update coefficients for the regular covariates
      if (iter==1){

        if(dist$family == "binomial"){
          GLM0 = glm_fit_stats(X, y, family = dist, intercept = FALSE, weights = wts,
                               control = glm.control(maxit = 5, epsilon = .Machine$double.eps))

        }else{
          GLM0 = glm_fit_stats(X, y, family = dist, intercept = FALSE, weights = wts)
        }

        beta0 <- GLM0$bb
        dev0 <- GLM0$dev
      }else{
        eta = Xj %*% as.vector(beta[[d]])

        if(dist$family == "binomial"){
          #GLM2 = glm_fit_stats_LA(cbind(X,eta), y, family = dist, intercept = F, weights = wts)

          # GLM2 = try(glmfit_priv(cbind(X,eta), y, dist = dist, const = "off", pwts = wts))
          # if (inherits(GLM1, "try-error")) break
          # if (inherits(GLM2, "try-error")) break
          #GLM2 = glmfit_priv(cbind(X,eta), y, dist = dist, const = "off", pwts = wts)
          GLM2 = glm_fit_stats(cbind(X,eta), y, family = dist, intercept = F, weights = wts,
                               control = glm.control(maxit = 5, epsilon = .Machine$double.eps))
        }else{
          GLM2 = glm_fit_stats(cbind(X,eta), y, family = dist, intercept = F, weights = wts)
        }
        betatmp <- GLM2$bb
        devtmp <- GLM2$dev
        glmstats <- GLM2$stats
        beta0 = betatmp[-length(betatmp)]

        # stopping rule
        diffdev = devtmp-dev0
        dev0 = devtmp
        if(abs(diffdev) < TolFun*(abs(dev0)+1)) break

        # update scale of array coefficients and standardize
        NormBeta = normalize(beta, rep(1, r)*betatmp[length(betatmp)])
        beta = NormBeta$beta
        beta.lambda = NormBeta$lambda

        for (jj in 1:d) {
          # bsxfun(@times,beta.U{j},(beta.lambda').^(1/d))
          for(ii in 1:r) beta[[jj]][, ii] = beta[[jj]][, ii]*(beta.lambda[ii]^(1/d))
        }
        beta.lambda = rep(1,r)
      }

      # cyclic update of the array coefficients
      eta0 = as.matrix(X) %*% beta0

      for (j in 1:d) {
        #if (j == 1) cumkr = matrix(rep(1, r), ncol = r)

        #create Xj
        if (exists("Md", mode = "list")){
          A1 = do.call(as.array, Md[j])
          loopindex = sort((1:d)[-j], decreasing = T)
          stop = length(loopindex)
          A21 = as.matrix(beta[[loopindex[1]]])
          t = 1
          while (t < stop) {
            t = t + 1
            A21 = khatri_rao(A21, beta[[loopindex[t]]])
          }

          #A2 = matrix(KTtoTensor(A21), ncol = 1)
          Xj = array(A1 %*% A21, dim = c(n, p[j]*r))
          #tcrossprod(A1, t(A21))

        }else{
          A1 = tenmat(TM,c(d+1,j))
          loopindex = sort((1:d)[-j], decreasing = T)
          stop = length(loopindex)
          A21 = as.matrix(beta[[loopindex[1]]])
          t = 1
          while (t < stop) {
            t = t + 1
            A21 = khatri_rao(A21, beta[[loopindex[t]]])
          }
          Xj = array(A1 %*% A21, dim = c(n, p[j]*r))
        }

        #run GLM
        if(dist$family == "binomial"){
          #GLM1 = glm_fit_stats_LA(cbind(Xj,eta0), y, family = dist, intercept = F, weights = wts)

          # GLM1 = try(glmfit_priv(cbind(Xj,eta0), y, dist = dist, const = "off", pwts = wts))
          # if (inherits(GLM1, "try-error")) break

          #GLM1 = glmfit_priv(cbind(Xj,eta0), y, dist = dist, const = "off", pwts = wts)
          GLM1 = glm_fit_stats(cbind(Xj,eta0), y, family = dist, intercept = F, weights = wts,
                               control = glm.control(maxit = 5, epsilon = .Machine$double.eps))
        }else{
          GLM1 = glm_fit_stats(cbind(Xj,eta0), y, family = dist, intercept = F, weights = wts)
        }
        #GLM1 = glmfit_priv(cbind(Xj,eta0), y, dist = dist, const = 'off', pwts = wts)
        betatmp <- GLM1$bb
        #GLM3 = glm.fit(cbind(Xj,eta0), y, family = dist, intercept = T, weights = wts)
        dummy <- GLM1$dev
        # glmstats[[j]] <- GLM1$stats
        #renew the parameter
        beta[[j]] = array(betatmp[-length(betatmp)], dim = c(p[j],r))
        eta0 = eta0*betatmp[length(betatmp)] ##t次的 coefficients for the regular
      }
    }


    # record if it has a smaller deviance
    # print(paste0("iter: ", iter))
    if (any(inherits(GLM1, "try-error"), inherits(GLM2, "try-error"))) next

    # print(paste0("Replicate: ", rep, " - ", dev0<dev_final))
    # print(paste0("dev_final: ", dev_final))
    if(dev0<dev_final){
      beta0_final = beta0
      beta_final = beta
      glmstats_final = glmstats
      dev_final = dev0
    }

    if(Display){
      print(paste0('replicate: ', rep))
      print(paste0('iterates: ', iter))
      print(paste0('deviance: ', dev0))
      print(paste0('beta0: ', beta0))
    }
  }

  # output BIC of the final model. Note deviance = -2*log-likelihood
  if(d == 2){
    glmstats_final$BIC = dev_final + log(n)*(r*(p[1]+p[2]-r)+p0)
  }else{
    glmstats_final$BIC = dev_final + log(n)*(r*(sum(p[1:(length(p)-1)])-d+1)+p0)
  }

  if(is.null(X)){
    mu = t(t(matrix(KTtoTensor(beta_final))) %*% matrix(M, ncol = n))
    beta0_final = NULL
  }else{
    mu = X %*% beta0_final + t(t(matrix(KTtoTensor(beta_final))) %*% matrix(M, ncol = n))
  }

  if(dist$family == "gaussian"){
    glmstats_final$yhat = mu
  }else{
    glmstats_final$yhat = 1/(1+exp(-mu))
  }

  return(list("beta0_final" = beta0_final,
              "beta_final" = beta_final,
              "glmstats_final" = glmstats_final,
              "dev_final" = dev_final))

}
