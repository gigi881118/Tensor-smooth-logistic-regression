TensorRegNet = function(X,M,y,r,dist,lambda,
                              #pentype,
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
                              wts  = NULL){
  # B0 = NULL; Display = "off"; BurninMaxIter = 20;
  # BurninTolFun = 1e-2; BurninReplicates = 5; PenaltyMaxIter = 50
  # PenaltyTolFun = 1e-4; warn = FALSE; wts  = NULL
  # B0 = B0reg

  ## setting parameter
  if(is.null(wts)) wts = rep(1, dim(X)[1])

  # check positivity of tuning parameter
  if(lambda==0) stop('lambda=0 (no penalization); call TensorReg instead')

  # check validity of rank r
  if(is.null(r)) r = 1

  # decide least squares or GLM model
  if(dist$family == "gaussian"){
    isglm = FALSE
  }else{
    isglm = TRUE
    # translate to model specifier for sparse regression
    if(dist$family == "binomial"){
      glmmodel = 'logistic'
    }else{
      if(dist$family == "poisson") glmmodel = 'loglinear'
    }
  }

  #----------- check dimensions ---------------
  if(is.null(X)){
    dimX = c(tail(dim(M), 1), 1)
    X = array(rep(1, prod(dimX)), dim = dimX)
  }

  # check dimensions
  n = dim(X)[1]
  p0 = dim(X)[2]
  d = length(dim(M)) - 1 #dimension of array variates
  p = dim(M) #sizes array variates
  if(p[length(p)] != n){
    stop("dimension of M does not match that of X!")
  }

  #convert M into a tensor T
  TM = M

  #------------- mode-d matricization ------------------
  #if space allowing, pre-compute mode-d matricization of TM
  mem_info <- gc()
  # total_mem <- as.numeric(memory.limit()) * 1024 * 1024 #turn megabytes to bytes
  # available_mem <- total_mem - mem_info[2]

  # CAUTION: may cause out of memory on Linux
  # if(d*(8*prod(dim(TM)))<.75*available_mem){
  #
  # }
  Md <- list()
  for (dd in 1:d) {
    Md[[dd]] = tenmat(TM, c(d+1, dd))
  }


  #--------- loose convergence criterion -------
  #Burn-in stage (loose convergence criterion)
  if(Display){
    cat(" ================ Burn-in stage ... ================")
  }

  # reduce tensor size for reliable estimation in burnin stage
  if(is.null(B0)){
    #no user-supplied start point
    if(dist$family == "gaussian"){
      shrink_factor = (n/5) / (r*sum(p[1:(length(p) - 1)]))
    }else{
      if(dist$family == "binomial"){
        shrink_factor = (n/10) / (r*sum(p[1:(length(p) - 1)]))
      }else{
        if(dist$family == "poisson") shrink_factor = (n/10) / (r*sum(p[1:(length(p) - 1)]))
      }
    }

    if(shrink_factor <= 1){
      KR = TensorReg(X,M,y,r,dist,
                       MaxIter = BurninMaxIter,
                       TolFun = BurninTolFun,
                       Replicates = BurninReplicates,
                       wts = wts)
      dummy = KR$beta0_final
      beta_burnin = KR$beta_final
    }else{
      targetdim = c(round(p[1:(length(p) - 1)]/shrink_factor), n)
      if(any(targetdim <= 1)) stop('Cannot make a reasonable B0,
      the data is not suitable for fitting the sparse regression.')

      M_reduce = interpolate(M, targetdim) ##array_resize: mode = 1

      #estimate at reduced dimension
      KR = TensorReg(X,M_reduce,y,r,dist,
                       MaxIter = BurninMaxIter,
                       TolFun = BurninTolFun,
                       Replicates = BurninReplicates,
                       wts = wts)
      dummy = KR$beta0_final
      beta_burnin = KR$beta_final
      #resize back to original dimension
      beta_burnin = interpolate(beta_burnin, p[-length(p)])

      # warm start from coarsened estimate
      KR = TensorReg(X,M,y,r,dist,
                       B0 = beta_burnin, ##B0 = beta_burnin
                       MaxIter = BurninMaxIter,
                       TolFun = BurninTolFun,
                       wts = wts)
      dummy = KR$beta0_final
      beta_burnin = KR$beta_final
    }
  }else{
    #user-supplied start point

    ##ndim of B
    p1 = c()
    if(class(B0) == "array"){
      p1 = dim(B0)
    }else{
      if(class(B0) == "list"){
        for (i in 1:length(B0)) {
          p1 = c(p1, nrow(B0[[i]]))
        }
      }
    }
    ndimB = length(p1)

    if(ndimB != d) stop('dimension of B0 does not match that of data!')

    # turn B0 into a tensor (if it's not)
    if(class(B0) == "numeric") B0 = as.tensor(B0)

    #resize to compatible dimension (if it's not)
    if(sum(p1 != p[-length(p)]) > 0){
      B0 = interpolate(B0, p[-length(p)]) ##edit
    }
    # perform CP decomposition if it's not a ktensor of correct rank
    if(class(B0) == "array" |
       (class(B0) == "list" & ncol(B0[[1]]) != r)) {
      B0tensor = cp(as.tensor(B0), num_components = r)
      B0 = B0tensor$U
    }
    beta_burnin = B0
  }

  if(Display) cat(' ========== Penalization stage ===========')

  glmstats = list()
  dev0 = Inf
  beta = beta_burnin
  for (iter in 1:PenaltyMaxIter) {
    # update regular covariate coefficients
    if (iter==1) {
      eta = tcrossprod(tenmat(TM,d+1), t(tenmat(KTtoTensor(beta),1:d)))
    }else{
      #eta = Xj %*% betatmp[-length(betatmp)]
      eta = tcrossprod(Xj, t(betatmp[-length(betatmp)]))
    }


    #run GLM
    if(dist$family == "binomial"){
      GLM1 = glm_fit_stats(cbind(X, eta), y, family = dist, intercept = F, weights = wts,
                           control = glm.control(maxit = 5, epsilon = .Machine$double.eps))
    }else{
      GLM1 = glm_fit_stats(cbind(X, eta), y, family = dist, intercept = F, weights = wts)
    }


    betatmp = GLM1$bb
    devtmp = GLM1$dev
    glmstats = GLM1$stats

    beta0 = betatmp[1:p0]
    # stopping rule
    diffdev = devtmp-dev0
    dev0 = devtmp

    # print(paste0("diff: ",  abs(diffdev), " penalty: ", PenaltyTolFun*(abs(dev0)+1)))
    if (abs(diffdev) < PenaltyTolFun*(abs(dev0)+1)) break

    if(!smooth){
      # update scale of array coefficients and standardize
      NormBeta = normalize(beta, rep(1,r)*betatmp[length(betatmp)])
      beta = NormBeta$beta
      beta.lambda =  NormBeta$lambda
      for (j in 1:d) {
        # bsxfun(@times,beta.U{j},(beta.lambda').^(1/d))
        for(i in 1:r) beta[[j]][, i] = beta[[j]][, i]*(beta.lambda^(1/d))[i]
      }
      beta.lambda = rep(1,r)
    }

    #cyclic update of array regression coefficients
    eta0 = tcrossprod(as.matrix(X), t(beta0))
    #---------- rewrite ----------
    for (j in 1:d) {
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
        Xj = array(tcrossprod(A1, t(A21)), dim = c(n, p[j]*r))

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
        Xj = array(tcrossprod(A1, t(A21)), dim = c(n, p[j]*r))
      }

      if(smooth){
        #---------------- LW　---------------
        TCV = matrix(1, nrow = d, ncol = r)
        for (jj in 1:d) {
          beta_jj = beta[[jj]]
          end_b = nrow(beta_jj)

          if(r != 1){
            TCV[jj, ] = apply(abs(beta_jj[2:end_b, ] - beta_jj[1:(end_b - 1), ]), 2, sum)
          }else{
            TCV[jj, ] = sum(abs(beta_jj[2:end_b, ] - beta_jj[1:(end_b - 1), ]))
          }
        }
        TCV[j, ] = NA
        LW = apply(TCV, 2, function(x) prod(x, na.rm = TRUE))
        # print(LW)

        if(any(round(LW, 10)==0)){
          print(LW)
          LW[which(round(LW, 10)==0)] = 10
          # stop("This is an error message")
        }
        Rj = matrix(0, p[j]-1, p[j])
        diag(Rj) = -1
        diag(Rj[, -1]) = 1
        Ar = diag(r)
        Dj = kronecker(LW*Ar, Rj) #LW*Ar
        D = cbind(Dj, 0)

        Arr = matrix(0, 1, ncol(D))
        Arr[, ncol(Arr)] = 1
        DD = rbind(D, Arr)

        L = lambda
        if(dist$family == "binomial"){


          shiv = ADMMadj5::admm.genlasso(A = cbind(Xj, eta0), b = y, D = DD,
                                         lambda = L, alpha = 1, rho = 1
                                          ,reltol = 0.01 , maxiter = 5
          )
          betatmp = shiv$x
          beta[[j]] = matrix(betatmp[-length(betatmp)],p[j],r)
          eta0 = eta0*betatmp[length(betatmp)]
          #-------------------------- stop ----------------------------------

        }else{

          shiv = ADMM::admm.genlasso(A = cbind(Xj, eta0), b = y, D = D,
                                     lambda = L, alpha = 1, maxiter = 5
          )
          betatmp = shiv$x
          beta[[j]] = matrix(betatmp[-length(betatmp)],p[j],r)
          eta0 = eta0*betatmp[length(betatmp)]
        }
      }else{
        # if(is.null(lambda)){
        #   cv_lsq = cv.glmnet(family = dist$family, cbind(Xj,eta0), y, weights = wts)
        #   L = cv_lsq$lambda.min
        # }else{
        #   L = lambda
        # }

        L = lambda

        lsq_sparse = glmnet(family = dist$family,
                            cbind(Xj,eta0), y, lambda = L, weights = wts,
                            penalty.factor = c(rep(1, p[j]*r), 0),
                            standardize = F, intercept = F, alpha = penparam
                            , control = glmnet.control(mxitnr = 5),
                            maxit = 1e+3
        )
        betatmp = lsq_sparse$beta
        beta[[j]] = matrix(betatmp[-length(betatmp)],p[j],r)
        eta0 = eta0*betatmp[length(betatmp)]
      }
    }

    if(Display){
      print(paste0('replicate: ', rep))
      print(paste0('iterates: ', iter))
      print(paste0('deviance: ', dev0))
      print(paste0('beta0: ', beta0))
    }
  }
  beta0_final = beta0
  beta_final = beta

  # find a scaling for the estimates
  beta_scale  <- lapply(1:d, function(j) matrix(0, ncol = r, nrow = p[j]))
  #eta0 = X %*% beta0
  eta0 = tcrossprod(as.matrix(X), t(beta0))
  for (j in 1:d) {
    idxj = 1:d
    idxj = idxj[-j]
    if (exists("Md", mode = "list")){
      A1 = do.call(as.array, Md[j])

      #A2
      idx = idxj[seq(length(idxj), 1, by = -1)]
      if(length(idx) == 1){
        A2 = beta[[idx]]
      }else{
        for (i in 1:length(idx)) {
          if(i == 1){A2 = beta[[idx[i]]]}else{
            #A2 = KhatriRao(A2, beta[[idx[i]]])
            A2 = khatri_rao(A2, beta[[idx[i]]])
          }
        }
      }

      #Xj = array(A1 %*% A2, dim = c(n, p[j]*r))
      Xj = array(tcrossprod(A1, t(A2)), dim = c(n, p[j]*r))
    }else{
      A1 = tenmat(TM,c(d+1,j))
      #A2
      idx = idxj[seq(length(idxj), 1, by = -1)]
      if(length(idx) == 1){
        A2 = beta[[idx]]
      }else{
        for (i in 1:length(idx)) {
          if(i == 1){A2 = beta[[idx[i]]]}else{
            #A2 = KhatriRao(A2, beta[[idx[i]]])
            A2 = khatri_rao(A2, beta[[idx[i]]])
          }
        }
      }

      #Xj = array(A1 %*% A2, dim = c(n, p[j]*r))
      Xj = array(tcrossprod(A1, t(A2)), dim = c(n, p[j]*r))
    }

    #GLM1 = glmfit_priv(cbind(Xj, eta0), y = y, dist = dist, const = "off")
    #GLM2 = glm_fit_stats(cbind(Xj, eta0), y = y, family = dist, intercept = FALSE)
    if(dist$family == "binomial"){
      GLM2 = glm_fit_stats(cbind(Xj, eta0), y = y, family = dist, intercept = FALSE,
                           control = glm.control(maxit = 5, epsilon = .Machine$double.eps))

      #GLM2 = glm_fit_stats_LA(cbind(Xj, eta0), y = y, family = dist, intercept = FALSE)
    }else{
      GLM2 = glm_fit_stats(cbind(Xj, eta0), y = y, family = dist, intercept = FALSE)
    }
    #GLM2 = glm_fit_stats(cbind(Xj, eta0), y = y, family = dist, intercept = FALSE)


    # glmstats[[d]] = GLM2$stats

    se = rep(0, p[j]*r+1)
    idx = which(!is.na(GLM2$bb))
    se[idx] = GLM2$stats$se
    #beta_scale{j} = reshape(glmstats{d}.se(1:end-1),p(j),r)
    beta_scale[[j]] = matrix(se[-length(se)], p[j], r)
  }

  #output the BIC
  cutoff = 1e-8
  # if(smooth){
  #   lowerBnd = log(.Machine$double.eps) ; upperBnd = -lowerBnd
  #   ilinkFun = function(eta) 1 / (1 + exp(-constrain(eta, lowerBnd,upperBnd)))
  #
  #   logit_piN = X%*%as.matrix(beta0_final) + apply(M, d+1, function(x)sum(x * KTtoTensor(beta_final)))
  #   piN <- ilinkFun(logit_piN)
  #   (dev0 = -2 * sum(y*log(piN) + (1-y)*log(1-piN)))
  # }

  if(d == 2){
    m1 = sum((abs(beta[[1]]) > cutoff)!= 0) + sum((abs(beta[[2]]) > cutoff)!= 0) - r*r
    glmstats$BIC = dev0 + log(n)*max(m1, 0)
  }else{
    if(smooth){
      # m1 = sum(unlist(lapply(beta, function(j){sum(unique(as.vector(j)) != 0)}))) - r*(d-1)
      m1 = sum(unlist(lapply(beta, function(j){sum(round(j, 4) != 0)}))) - r*(d-1)
    }else{
      # m1 = sum(unlist(lapply(beta, function(j){sum(j != 0)}))) - r*(d-1)
      m1 = sum(unlist(lapply(beta, function(j){sum(round(j, 4) != 0)}))) - r*(d-1)
    }
    glmstats$BIC = dev0 + log(n)*max(m1, 0)

  }

  if(is.null(X)){
    mu = t(t(matrix(KTtoTensor(beta_final))) %*% matrix(M, ncol = n))
    beta0_final = NULL
  }else{
    mu = X %*% beta0_final + t(t(matrix(KTtoTensor(beta_final))) %*% matrix(M, ncol = n))
  }

  if(dist$family == "gaussian"){
    glmstats$yhat = mu
  }else{
    glmstats$yhat = 1/(1+exp(-mu))
  }

  # say goodbye
  if(Display){
    print(" DONE! ")
  }

  return(list("beta0_final" = beta0_final,
              "beta_final" = beta_final,
              "beta_scale" = beta_scale,
              "glmstats" = glmstats))
}
