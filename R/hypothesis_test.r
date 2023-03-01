


# library('CompQuadForm')
# library('dplyr')


# 
# Carter, Schnepel, and Steigerwald test
# 

#' Eigenvalues for CSS Test
#'
#' This function gets the eigenvalues for a test based on
#' Carter, Schnepel, and Steigerwald (2017).
#'
#' @param xg List of cluster-specific submatrices of design matrix
#' @param xtx Inverse of design matrix squared
#' @param c.0 Linear hypothesis vector
#' @return A vector of eigenvalues
css.test.lambdas <- function(xg, xtx, c.0){
  # primitives
  MU = xtx %*% c.0
  RHO_G = Map(function(x) (-1 / (nrow(x))), x=xg)
  ETA_G = Map(function(x, rho) ((t(x) %*% x)), x=xg, rho=RHO_G) # eta
  LAM_0 = Map(function(eta, rho) ((t(MU) %*% eta %*% MU)), eta=ETA_G, rho=RHO_G) # eta
  return(unlist(LAM_0))
}


# 
# Bell and McCaffrey test
# 

#' Eigenvalues for BM Test
#'
#' This function gets the eigenvalues for a test based on
#' Bell and McCaffrey (2002).
#'
#' @param xg xg List of cluster-specific submatrices of design matrix
#' @param xtx Inverse of design matrix squared
#' @param ag List of adjustment matrices
#' @param c.0 Linear hypothesis vector
#' @return A vector of eigenvalues
bm.test.lambdas <- function(xg, xtx, ag, c.0){
  # primitives
  MU = xtx %*% c.0
  DELTA_G = Map(function(x, a) (t(x) %*% t(a) %*% x), x=xg, a=ag) # delta
  ETA_G = Map(function(x, rho) ((t(x) %*% x)), x=xg, rho=RHO_G) # eta
  
  MUPHIMU_G = Map(function(x, a){
    mpm.0 = (t(a) %*% (x %*% MU))
    mpm.1 = t(mpm.0) %*% mpm.0
    return(mpm.1)
  }, x=xg, a=ag) # mu phi mu
  
  # transformations
  V.0 = xtx %*% Reduce('+', ETA_G) %*% xtx # v
  Delta.0 = Map(function(delta) (delta %*% MU), delta=DELTA_G) # delta mu
  Theta.0 = Map(function(delta) (xtx %*% (delta) %*% MU),
                delta=DELTA_G) # XTX theta mu
  
  # big matrices
  Delta.1 = do.call(cbind, Delta.0) # Delta
  Theta.1 = do.call(cbind, Theta.0) # Omega
  Phi.1 = diag(unlist(MUPHIMU_G))
  
  DOD.0 = (t(Delta.1) %*% V.0 %*% Delta.1) + Phi.1 - ((t(Delta.1) %*% Theta.1) + (t(Theta.1) %*% Delta.1))
  DOD.1 = zapsmall(DOD.0)
  lam.0 = eigen(DOD.0)$values
  
  if(sum(abs(Im(lam.0))) > (max(abs(Re(lam.0))) * 10^(-10))){
    stop('bm.test.lambdas(): complex eigenvalues error')
  }else{
    lam.0 = Re(lam.0)
  }
  
  return(lam.0)
}


# 
# Meiselman test
# 

#' Eigenvalues for Meiselman Test
#'
#' This function gets an intermediate thing for the eigenvalues
#' for a test based on Meiselman (2021).
#'
#' @param xg List of cluster-specific submatrices of design matrix
#' @param xtx Inverse of design matrix squared
#' @param c.0 Linear hypothesis vector
#' @return A list containing two principal submatrices of the key matrix
meis.test.lambdas.1 <- function(xg, xtx, ag, c.0){
  # primitives
  MU = xtx %*% c.0
  DELTA_G = Map(function(x, a) (t(x) %*% t(a) %*% x), x=xg, a=ag) # delta
  MUPHIMU_G = Map(function(x, a){
    mpm.0 = (t(a) %*% (x %*% MU))
    mpm.1 = t(mpm.0) %*% mpm.0
    return(mpm.1)
  }, x=xg, a=ag) # mu phi mu
  
  # transformations
  Delta.0 = Map(function(delta) (t(MU) %*% delta), delta=DELTA_G)
  
  # big matrices
  Delta.1 = do.call(rbind, Delta.0) # Delta
  Phi.1 = diag(unlist(MUPHIMU_G))
  DOD.0 = Phi.1 - (Delta.1 %*% xtx %*% t(Delta.1))
  
  ret = list()
  ret$dod = DOD.0
  ret$v = t(c.0) %*% xtx %*% c.0
  return(ret)
}

#' Eigenvalues for Meiselman Test
#'
#' This function gets the eigenvalues for a test based on
#' Meiselman (2021).
#'
#' @param m1 A list containing two principal submatrices of the key matrix
#' @param c.0 Linear hypothesis vector
#' @param q Desired quantile
#' @return A vector of eigenvalues
meis.test.lambdas.2 <- function(m1, c.0, q){
  w.13 = -(m1$v) / (q ^ 2) # -fT Om f / (z^2)
  lam.0 = c(eigen(m1$dod)$values, w.13)
  return(lam.0)
}

#' P-Value for Meiselman Test
#'
#' This function gets the p-value for a test based on
#' Meiselman (2021).
#'
#' @param m1 A list containing two principal submatrices of the key matrix
#' @param c.0 Linear hypothesis vector
#' @param q Desired quantile
#' @return Probability that abs(test statistic) is greater than q
p.meis <- function(m1, c.0, q){
  lam.1 = meis.test.lambdas.2(m1, c.0, q)
  if(sum(abs(Im(lam.1))) > (max(abs(Re(lam.1))) * 10^(-10))){
    stop('p.meis(): complex eigenvalues error')
  }else{
    lam.1 = Re(lam.1)
  }
  lam.2 = lam.1 / sum(abs(lam.1))
  p.ret = 1 - CompQuadForm::imhof(0, lam.2)$Qq
  return(p.ret)
}

#' Critical Value for Meiselman Test
#'
#' This function gives the critical value for a
#' test based on Meiselman (2021).
#'
#' @param m1 A list containing two principal submatrices of the key matrix
#' @param c.0 Linear hypothesis vector
#' @param q Desired quantile
#' @return Critical value for a test with size alpha
q.meis <- function(m1, c.0, alpha=0.05){
  f1 <- function(q){
    p.curr = p.meis(m1, c.0, q)
    p.diff = p.curr - alpha
    return(p.diff)
  }
  q.ret = stats::uniroot(f1, c(0.01, 10), tol=(10^(-6)), extendInt='downX')$root
  
  return(q.ret)
}



# slower but more straightforward, for testing

#' Eigenvalues for Meiselman Test
#'
#' This function gets the eigenvalues for a test based on
#' Meiselman (2021).
#'
#' @param xg List of design matrix, separated by cluster
#' @param xtx Inverse of design matrix squared
#' @param c.0 Linear hypothesis vector
#' @return A vector of eigenvalues
meis.test.lambdas.slow <- function(xtx, ag, c.0, t.0, omg, xs, ihg, ind.hi, ind.lo,
                                   D0=NULL){
  X_g = lapply(seq_len(length(ind.hi)), function(x) as.matrix(xs[ind.lo[x]:ind.hi[x],]))

  f.0 = xs %*% (xtx %*% c.0)
  f.1 = f.0 / t.0
  
  if(is.null(D0)){
    D_0 <- Map(function(x, a, ih) (t(ih) %*% (t(a) %*% (x %*% (xtx %*% c.0)))),
               x=X_g, a=ag, ih=ihg)
    D_1 = do.call(cbind, D_0)
  }else{
    D_1 = D0
  }
  D_2 = cbind(f.1, D_1)
  D_3 = cbind(-f.1, D_1)

  D_4 = lapply(seq_len(length(ind.hi)), function(x) as.matrix(D_2[ind.lo[x]:ind.hi[x],]))
  D_5 = lapply(seq_len(length(ind.hi)), function(x) as.matrix(D_3[ind.lo[x]:ind.hi[x],]))
  D_6 = Map(function(d4, d5, om) (t(d4) %*% om %*% d5),
            d4=D_4, d5=D_5, om=omg)
  D_7 = Reduce('+', D_6)
  ll = eigen(D_7)$values
  return(ll)
}

#' Test Parameters
#'
#' This function tests the validity of several basic
#' model parameters.
#'
#' @param data Data frame or matrix containing all outcomes, covariates, and group ids
#' @param x.names Vector of covariate names
#' @param cluster.id Name of cluster id variable
#' @param weights Name of weight variable, if applicable
#' @param c0 Linear hypothesis vector, either c0 or test.var must be supplied
#' @param test.var Index of variable to test, either c0 or test.var must be supplied
#' @param fe.id Name of fixed effect group id variable (defaults to cluster.id)
#' @param calling.function Function whose parameters are being checked
#' @return TRUE if no errors, otherwise does not return
test.params <- function(data, x.names, cluster.id, weights=NULL, c0=NULL,
                        test.var=NULL, fe.id=cluster.id, calling.function='unknown'){
  #
  # testing all inputs
  #
  
  # does data have covariates
  if(!all(x.names %in% colnames(data))){
    stop(paste(calling.function, 'error: x not all in data'))
  }
  # does data have cluster ids
  if(!(cluster.id %in% colnames(data))){
    stop(paste(calling.function, 'error: cluster.id not in data'))
  }
  # does data have fixed effect group ids
  if(!(fe.id %in% colnames(data))){
    stop(paste(calling.function, 'error: fe.id not in data'))
  }
  # does data have weights
  if(!is.null(weights)){
    if(!(weights %in% colnames(data))){
      stop(paste(calling.function, 'error: weights not in data'))
    }
  }
  
  
  # is fe.id nested within cluster.id
  if(fe.id != cluster.id){
    tmp1 = length(unique(data[, fe.id]))
    tmp2 = paste('tag', data[, cluster.id], data[, fe.id], 'end', sep='.')
    tmp3 = length(unique(tmp2))
    if(tmp1 != tmp3){
      stop(paste(calling.function, "error: fe.id is not nested inside cluster.id"))
    }
  }
  
  # is either test.var or c0 filled in, and is one correct
  if(!is.null(c0)){
    if(!is.numeric(c0)){
      stop(paste(calling.function, "error: c0 must be numeric"))
    }
    if(length(c(c0)) != length(x.names)){
      stop(paste(calling.function, "error: c0 must be K-vector"))
    }
    if(sum(abs(c0)) == 0){
      stop(paste(calling.function, "error: c0 cannot be 0"))
    }
  }else{
    if(!is.null(test.var)){
      if(!(test.var %in% x.names)){
        stop(paste(calling.function, "error: test.var not in x.names"))
      }
      # if(!is.numeric(test.var)){
      #   stop(paste(calling.function, "error: test.var must be numeric"))
      # }
      # if((test.var < 1) | (test.var > length(x.names))){
      #   stop(paste(calling.function, "error: test.var must refer to one of K covariates"))
      # }
    }else{
      stop(paste(calling.function, "error: c0 or test.var must be included"))
    }
  }
  
  if(!(crve.type %in% c('CR0', 'CR1', 'CR2', 'CR3'))){
    stop("p.value.meis() error: crve.type should be in (CR0, CR1, CR2, CR3)")
  }
  
  return(T)
}

#' Absorb Fixed Effects and Reorder Data
#'
#' This function reorganizes the data by cluster and fixed effect,
#' absorbing group fixed effects and applying weights.
#'
#' @param data Data frame or matrix containing all outcomes, covariates, and group ids
#' @param x.names Vector of covariate names
#' @param cluster.id Name of cluster id variable
#' @param weights Name of weight variable, if applicable
#' @param fe.id Name of fixed effect group id variable (defaults to cluster.id)
#' @param crve Type of cluster-robust variance estimator (CR0, CR1, CR2, CR3)
#' @return A list containing: Xdw.h (covariates by cluster),
#'         XTX.inv (inverse squared design matrix),
#'         A.h (adjustment matrices by cluster),
#'         Ydw.h (outcomes by cluster, if applicable)
absorb.reorg <- function(data, x.names, cluster.id, weights, fe.id, crve, y.name=NULL){
  #
  # reorder design matrix
  #
  data.0 = data
  if(!is.null(weights)){
    data.0$weights = data.0[, weights]
  }else{
    data.0$weights = 1
  }
  data.0$clu.var = data.0[, cluster.id]
  data.0$fe.var = data.0[, fe.id]
  data.0 = data.0[with(data.0, order(clu.var, fe.var)), ]
  data.0$id = 1:nrow(data.0)
  
  X = as.matrix(data.0[, x.names, drop=F])
  W = as.matrix(data.0[, 'weights', drop=F])
  
  #
  # Absorb fixed effects
  #
  
  # groups for fixed effects
  data.0$one = 1
  groups.0 = data.0 %>%
    group_by(clu.var, fe.var) %>%
    summarize(obs=sum(one),
              id.min=min(id),
              id.max=max(id))
  groups.0 = as.data.frame(groups.0)
  
  # groups for clustering
  groups.1 = data.0 %>%
    group_by(clu.var) %>%
    summarize(obs=sum(one),
              id.min=min(id),
              id.max=max(id))
  groups.1 = as.data.frame(groups.1)
  
  ind.lo = groups.0$id.min
  ind.hi = groups.0$id.max
  ngrp = nrow(groups.0)
  
  X.g = lapply(seq_len(ngrp), function(x) as.matrix(X[ind.lo[x]:ind.hi[x],, drop=F]))
  W.g = lapply(seq_len(ngrp), function(x) as.matrix(W[ind.lo[x]:ind.hi[x],, drop=F]))
  
  # absorb FE and apply weights
  Xdw.g = Map(function(x, w){
    x.new = x
    for(j in 1:ncol(x)){
      x.new[, j] = sqrt(w) * (x[, j] - (sum(w * x[, j]) / sum(w)))
    }
    return(x.new)
  }, x=X.g, w=W.g)
  Xdw = do.call(rbind, Xdw.g)
  
  
  #
  # Generate parameters for meis.test.lambdas.1()
  #
  XTX = t(Xdw) %*% Xdw
  XTX.inv = solve(XTX)
  
  ind.lo = groups.1$id.min
  ind.hi = groups.1$id.max
  nclu = nrow(groups.1)
  Xdw.h = lapply(seq_len(nclu), function(x) as.matrix(Xdw[ind.lo[x]:ind.hi[x],, drop=F]))
  
  if(crve == 'CR0'){
    A.h = Map(function(x) (diag(nrow(x))), x=Xdw.h)
  }else if(crve == 'CR2'){
    A.h = adjustment.matrices(Xdw.h, XTX.inv, delta=0.5)
  }else if(crve == 'CR3'){
    A.h = adjustment.matrices(Xdw.h, XTX.inv, delta=1)
  }else{
    stop("p.value.meis() error: wrong crve.type")
  }
  
  model.0 = list()
  model.0$Xdw.h = Xdw.h
  model.0$XTX.inv = XTX.inv
  model.0$A.h = A.h
  
  if(!is.null(y.name)){
    Y = as.matrix(data.0[, y.name, drop=F])
    
    ind.lo = groups.0$id.min
    ind.hi = groups.0$id.max
    ngrp = nrow(groups.0)
    
    Y.g = lapply(seq_len(ngrp), function(x) as.matrix(Y[ind.lo[x]:ind.hi[x],, drop=F]))
    
    # absorb FE and apply weights
    Ydw.g = Map(function(x, w){
      x.new = x
      for(j in 1:ncol(x)){
        x.new[, j] = sqrt(w) * (x[, j] - (sum(w * x[, j]) / sum(w)))
      }
      return(x.new)
    }, x=Y.g, w=W.g)
    Ydw = do.call(rbind, Ydw.g)
    
    ind.lo = groups.1$id.min
    ind.hi = groups.1$id.max
    nclu = nrow(groups.1)
    Ydw.h = lapply(seq_len(nclu), function(x) as.matrix(Ydw[ind.lo[x]:ind.hi[x],, drop=F]))
    
    model.0$Ydw.h = Ydw.h
  }
  
  return(model.0)
}

#' Calculate P-Value
#'
#' This function calculates a p-value for the test statistic t0
#' according to Meiselman (2021).
#'
#' @param t0 The test statistic
#' @param data Data frame or matrix containing all outcomes, covariates, and group ids
#' @param x.names Vector of covariate names
#' @param cluster.id Name of cluster id variable
#' @param weights Name of weight variable, if applicable
#' @param c0 Linear hypothesis vector, either c0 or test.var must be supplied
#' @param test.var Index of variable to test, either c0 or test.var must be supplied
#' @param fe.id Name of fixed effect group id variable (defaults to cluster.id)
#' @param crve.type Type of cluster-robust variance estimator (CR0, CR1, CR2, CR3)
#' @return numeric, the probability that a test statistic with a magnitude of abs(t0) or greater is observed
#' @export
p.value.meis <- function(t0, data, x.names, cluster.id, weights=NULL, c0=NULL, test.var=NULL,
                         fe.id=cluster.id, crve.type='CR1'){
  #
  # testing all inputs
  #
  test.params(data, x.names, cluster.id, weights=weights, c0=c0,
              test.var=test.var, fe.id=fe.id, calling.function='p.value.meis()')
  
  if(!is.null(c0)){
    c1 = c(c0)
  }else if(!is.null(test.var)){
    c1 = rep(0, length(x.names))
    c1[x.names == test.var] = 1
  }
  
  if(crve.type %in% c('CR0', 'CR2', 'CR3')){
    t1 = t0
    crve = crve.type
  }else if(crve.type == 'CR1'){
    nclu = length(unique(cluster.id))
    t1 = t0 * sqrt(nclu / (nclu - 1))
    crve = 'CR0'
  }
  
  model.0 = absorb.reorg(data, x.names, cluster.id, weights, fe.id, crve)
  
  #
  # pass inputs to other functions
  #
  c2 = t(t(c1))
  m1 = meis.test.lambdas.1(model.0$Xdw.h, model.0$XTX.inv, model.0$A.h, c2)
  p0 = p.meis(m1, c2, t1)
  return(p0)
}


#' Calculate Confidence Interval
#'
#' This function calculates a confidence interval for the
#' coefficient (or linear combination of coefficients) defined
#' by test.var (or c0) according to Meiselman (2021).
#'
#' @param data Data frame or matrix containing all outcomes, covariates, and group ids
#' @param y.name Name of dependent variable
#' @param x.names Vector of covariate names
#' @param cluster.id Name of cluster id variable
#' @param alpha Size of test associated with the confidence interval (e.g. a test of size 0.05 is associated with a 95% confidence interval)
#' @param weights Name of weight variable, if applicable
#' @param c0 Linear hypothesis vector, either c0 or test.var must be supplied
#' @param test.var Index of variable to test, either c0 or test.var must be supplied
#' @param fe.id Name of fixed effect group id variable (defaults to cluster.id)
#' @param crve.type Type of cluster-robust variance estimator (CR0, CR1, CR2, CR3)
#' @return numeric, the lower and upper bound of the confidence interval
#' @export
conf.interval.meis <- function(data, y.name, x.names, cluster.id, alpha=0.05,
                         weights=NULL, c0=NULL, test.var=NULL,
                         fe.id=cluster.id, crve.type='CR1'){
  #
  # testing all inputs
  #
  test.params(data, x.names, cluster.id, weights=weights, c0=c0,
              test.var=test.var, fe.id=fe.id, calling.function='p.value.meis()')
  
  # does data have dependent variable
  if(!(y.name %in% colnames(data))){
    stop(paste(calling.function, 'error: y.anme not in data'))
  }
  
  if(!is.null(c0)){
    c1 = c(c0)
  }else if(!is.null(test.var)){
    c1 = rep(0, length(x.names))
    c1[x.names == test.var] = 1
  }
  
  crve = crve.type
  if(crve.type == 'CR1'){
    crve = 'CR0'
  }
  
  if(!is.numeric(alpha)){
    stop("conf.interval.meis() error: alpha must be a number")
  }
  if((alpha <= 0) | (alpha >= 1)){
    stop("conf.interval.meis() error: alpha must be in (0,1)")
  }
  
  model.0 = absorb.reorg(data, x.names, cluster.id, weights, fe.id, crve)
  Xdw = do.call(rbind, model.0$Xdw.h)
  Ydw = do.call(rbind, model.0$Ydw.h)
  # OLS
  XTY = t(Xdw) %*% Ydw
  B.hat = model.0$XTX.inv %*% XTY
  # CRVE
  Eh_h = Map(function(x, y) (y - x %*% B.hat), x=Xdw.h, y=Ydw.h)
  VAR_0 = Map(function(x, eh) ((t(x) %*% eh) %*% (t(eh) %*% x)), x=Xdw.h, eh=Eh.h)
  VAR_1 = Reduce('+', VAR_0)
  VAR_2 = model.0$XTX.inv %*% VAR_1 %*% model.0$XTX.inv
  
  #
  # pass inputs to other functions
  #
  c2 = t(t(c1))
  m1 = meis.test.lambdas.1(model.0$Xdw.h, model.0$XTX.inv, model.0$A.h, c2)
  q.star = q.meis(m1, c2, alpha=alpha)
  
  #
  # Construct confidence interval
  #
  lb = (t(c2) %*% B.hat) - (q.star * (t(c2) %*% VAR_2 %*% c2))
  ub = (t(c2) %*% B.hat) + (q.star * (t(c2) %*% VAR_2 %*% c2))
  ci.0 = c(lower.bound=lb, upper.bound=ub)
  
  return(ci.0)
}

