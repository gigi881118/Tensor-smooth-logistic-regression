### auxiliary functions

#------- constrain ---------
#Constrain between upper and lower limits, and do not ignore NaN
constrain <- function(x,lower,upper){
x[x<lower] = lower
x[x>upper] = upper
return(x)
}


#---- glm fit stats ---------
#glm regression used for TR reg.
glm_fit_stats = function(x, y, weights = rep(1, nobs), start = NULL,
                         etastart = NULL, mustart = NULL, offset = rep(0, nobs),
                         family = gaussian(), control = list(), intercept = FALSE,
                         estdisp = 'on'){

  nobs <- NROW(y)

  ## define stats needed for glmpriv
  N = NULL #needed only for binomial
  y2 = y
  dist = family
  if(dist$family == "binomial"){
    if(dim(y2)[2] == 1){
      # N will get set to 1 below
      if(any(y2 < 0 || y > 1)){
        stop("stats:glmfit:BadData,
           For the binomial distribution, Y must be a binary vector or
           \na matrix with two columns with the number of trials in the second column.")
      }
    }else{
      if(dim(y2)[2] == 2){
        y2[which(y2[, 2] == 0), 2] = NaN
        N = y2[, 2]
        y2 = y2[, 1] / N
        if(any(y2 < 0 || y2 > 1)){
          stop("stats:glmfit:BadData,
           For the binomial distribution, Y must be a binary vector or
           \na matrix with two columns with the number of trials in the second column.")
        }
      }else{
        stop('stats:glmfit:MatrixOrBernoulliRequired,
           Y must be a two column matrix or a vector for the binomial distribution.')
      }
    }
  }
  pwts = weights
  #Remove missing values from the data.  Also turns row vectors into columns.
  rmNaN = statRmNan(list(y, x, offset, pwts, N))
  anybad = rmNaN$anybad
  wasnan = rmNaN$wasnan
  switch (anybad,
          "2" = stop(paste0('stats:glmfit:InputSizeMismatch', 'Number of observations in X and Y must match.')),
          "3" = stop(paste0('stats:glmfit:InputSizeMismatch', 'Lengths of OFFSET and Y must match.')),
          "4" = stop(paste0('stats:glmfit:InputSizeMismatch', 'Lengths of PWTS and Y must match.'))
  )

  ## define weights and offset if needed
  #if (is.null(weights)) weights <- rep.int(1, nobs)
  #if (is.null(offset))  offset <- rep.int(0, nobs)

  ##pwts, offset, N
  if(is.null(pwts)){
    pwts = 1
  }else{
    # A zero weight means ignore the observation, so n is reduced by one.
    # Residuals will be computed, however.
    if(any(pwts == 0)) n = n - sum(pwts == 0)
  }
  if(is.null(offset)) offset = 0
  if(is.null(N)) N = 1

  #print(paste0("intercept: ", intercept))

  #----------- fit GLM -----------------------------------------
  GLM = glm.fit(x = x, y = y, weights = pwts, start = start,
                etastart = etastart, mustart = mustart, offset = offset,
                family = family, control = control, intercept = intercept)

  stats = GLM
  names(stats)[names(stats)=="coefficients"] <- "beta"
  stats$fitted.values <- NULL
  stats$deviance <- NULL

  bb = GLM$coefficients
  bb[is.na(bb)] = 0
  dev = GLM$deviance

  # stats = list()
  stats$beta = bb

  # dfe = GLM$df.residual
  # stats$dfe = dfe
  # mu = GLM$fitted.values
  # y1 = GLM$y
  # R = GLM$R

  if(intercept == T){
    p = dim(x)[2] + 1
  }else{
    p = dim(x)[2]
  }

  perm = GLM$qr$pivot
  perm2 = perm[1:GLM$rank]
  s1 = summary.glm(GLM)

  stats$se = rep(0, p)
  stats$se[perm2] = s1$coefficients[, 2]

  return(invisible(list("bb" = bb,
                        "dev" = dev,
                        "stats" = stats)))


}

#------ insertnan ---------
#insertnan
insertnan <- function(wasnan, x){
  defaultW <- getOption("warn")
  options(warn = -1)
  if(!any(wasnan)) return(x); options(warn = defaultW)

  ok=!wasnan
  len = length(wasnan)
  if(!any(class(x) == "matrix")) x = as.matrix(x)
  y1 = x
  if((dim(x)[1] == 1) & (sum(ok) > 1)) y1 = t(y1)
  p = tail(dim(y1), 1)
  if(is.character(y1) || is.list(y1)){
    x = matrix(" ", len, p)
  }else{
    if(is.numeric(y1)){x = matrix(NaN, nrow = len, ncol = p)}
  }
  x[ok, ] = y1; options(warn = defaultW)
  return(x)

}

#---------- normalize ---------
normalize <- function(x, lambda){
  for (r in 1:length(lambda)) {
    for (n in 1:length(x)) {
      tmp = norm(as.matrix(x[[n]][, r]), type = "F")
      if (tmp > 0) x[[n]][, r] = x[[n]][, r] / tmp
      lambda[r] = lambda[r] * tmp
    }
  }
  # Check that all the lambda values are positive
  idx = which(lambda < 0)
  x[[1]][, idx] = -1 * x[[1]][, idx]
  lambda[idx] = -1 * lambda[idx]


  # Sort
  if(length(table(sort(lambda))) == 1){
    return(list("beta" = x, "lambda" = lambda))
  }else{
    idx2 = order(lambda, decreasing = T)
    a = list()
    for (r in 1:length(x)) a[[r]] <- x[[r]][, idx2]

    return(list("beta" = a, "lambda" = lambda[idx2]))
  }
}


#-------- statRmNan ---------
statRmNan <- function(x){
  idx1 = unlist(lapply(x, is.null))
  notNull = which(idx1 == FALSE)
  data = x[notNull]
  if(length(unique(unlist(lapply(x, nrow)))) != 1) stop("Error: row of data is not equal")
  if(any(class(data[[1]]) == "numeric")){n = length(data[[1]])}else{n = nrow(data[[1]])}
  X1 = matrix(unlist(data), nrow = n)
  #anybad
  anybad = c()
  for (i in notNull) {
    if(any(is.na(x[[i]]))) anybad = c(anybad, i)
  }
  if(length(anybad) == 0) anybad = 0
  #wasnan
  wasnan = apply(X1, 1, function(x) any(is.na(x)))
  if(sum(wasnan) == 0){
    return(list("anybad" = anybad,
                "wasnan" = wasnan,
                "var" = x))
  }else{
    idx = which(wasnan == TRUE)
    for (i in notNull) {
      if(any(class(x[[i]]) == "numeric")){
        x[[i]] = x[[i]][-idx]
      }else{
        x[[i]] = x[[i]][-idx, ]
      }
    }
    return(list("anybad" = anybad,
                "wasnan" = wasnan,
                "var" = x))
  }
}

#--------- tenmat ----------
tenmat <- function(TM, rdims, cdims = NULL){
  if (max(rdims) > length(dim(TM))) {
    stop("dimension of M does not match that of rdims", call. = FALSE)
  }
  #define index
  cdims = (1:length(dim(TM)))[-rdims]
  c_Dim = prod(dim(TM)[cdims])
  r_Dim = prod(dim(TM)[rdims])
  #dims = c(rdims, cdims)
  #D = length(dim(TM))

  outputM = array(aperm(TM, c(rdims, cdims)), dim = c(r_Dim, c_Dim))
  return(outputM)
}


