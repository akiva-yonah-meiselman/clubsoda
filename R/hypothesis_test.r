


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
css.test.gammas <- function(xg, xtx, c.0){
  MU = xtx %*% c.0
  gamma = Map(function(x){
    return(t(MU) %*% (t(x) %*% x) %*% MU)
  }, x=xg)
  
  return(unlist(gamma))
}

#' Degrees of freedom for CSS Test
#'
#' This function gets G*, the degrees of freedom for a test based on
#' Carter, Schnepel, and Steigerwald (2017).
#'
#' @param xg List of cluster-specific submatrices of design matrix
#' @param xtx Inverse of design matrix squared
#' @param c.0 Linear hypothesis vector
#' @return A scalar degrees of freedom
css.test.gstar <- function(xg, xtx, c.0){
  gamma = css.test.gammas(xg, xtx, c.0)
  gstar = (sum(gamma)^2) / sum(gamma^2)
  return(gstar)
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
#' @param ag List of adjustment matrices, either from Bell and McCaffrey (2002) or from Niccodemi et al. (2020)
#' @param c.0 Linear hypothesis vector
#' @return A vector of eigenvalues
bm.test.lambdas <- function(xg, xtx, ag, c.0, atype='NAAMW', xgtxg=NULL){
  MU = xtx %*% c.0
  
  if(atype == 'NAAMW'){
    if(is.null(xgtxg)){
      xgtxg = Map(function(x) (t(x) %*% x), x=xg) # eta
    }
    delta_g = Map(function(x, a) (x %*% a %*% MU), x=xgtxg, a=ag)
    omega_g = Map(function(d, a) (t(MU) %*% t(a) %*% d), d=delta_g, a=ag)
  }else if(atype == 'BM'){
    # xgtaaxg = Map(function(x, a) (t(x) %*% a %*% a %*% x), x=xg, a=ag)
    delta_g = Map(function(x, a) (t(x) %*% a %*% x %*% MU), x=xg, a=ag)
    omega_g = Map(function(x, a) (t(MU) %*% t(x) %*% a %*% a %*% x %*% MU), x=xg, a=ag)
  }else{
    stop('need proper atype')
  }
  
  Delta.0 = do.call(cbind, delta_g)
  Omega.0 = diag(unlist(omega_g))
  
  ODD.0 = Omega.0 - (t(Delta.0) %*% xtx %*% Delta.0)
  
  lambda.0 = eigen(ODD.0)$values
  
  if(sum(abs(Im(lambda.0))) > (max(abs(Re(lambda.0))) * 10^(-10))){
    stop('bm.test.lambdas(): complex eigenvalues error')
  }else{
    lambda.0 = Re(lambda.0)
  }
  
  return(lambda.0)
}



#' Degrees of freedom for BM Test
#'
#' This function gets m, the degrees of freedom for a test based on
#' Bell and McCaffrey (2002).
#'
#' @param xg xg List of cluster-specific submatrices of design matrix
#' @param xtx Inverse of design matrix squared
#' @param ag List of adjustment matrices, either from Bell and McCaffrey (2002) or from Niccodemi et al. (2020)
#' @param c.0 Linear hypothesis vector
#' @return A vector of eigenvalues
bm.test.df <- function(xg, xtx, ag, c.0, atype='BM', xgtxg=NULL){
  lam = bm.test.lambdas(xg, xtx, ag, c.0, atype=atype, xgtxg=xgtxg)
  m = (sum(lam)^2) / sum(lam^2)
  return(m)
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
p.meis <- function(xtx, c.0, q,
                   lambdas=NULL, design.info=NULL){
  if(is.null(lambdas)){
    if(is.null(design.info)){
      stop('need lambdas or design.info')
    }else{
      lambdas = bm.test.lambdas(design.info$xg, xtx, design.info$ag, c.0,
                                atype=design.info$atype, xgtxg=design.info$xgtxg)
    }
  }
  
  lam.0 = -t(c.0) %*% xtx %*% c.0 / (q^2)
  lam.1 = c(lam.0, lambdas)
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
q.meis <- function(xtx, c.0, alpha=0.05,
                   lambdas=NULL, design.info=NULL){
  f1 <- function(q){
    p.curr = p.meis(xtx, c.0, q,
                    lambdas=lambdas, design.info=design.info)
    p.diff = p.curr - alpha
    return(p.diff)
  }
  q.ret = stats::uniroot(f1, c(0.01, 10), tol=(10^(-6)), extendInt='downX')$root
  
  return(q.ret)
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
  clu.grps = data.0 %>%
    group_by(clu.var) %>%
    summarize(obs=sum(one),
              id.min=min(id),
              id.max=max(id))
  clu.grps = as.data.frame(clu.grps)
  
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
  
  ind.lo = clu.grps$id.min
  ind.hi = clu.grps$id.max
  nclu = nrow(clu.grps)
  Xdw.h = lapply(seq_len(nclu), function(x) as.matrix(Xdw[ind.lo[x]:ind.hi[x],, drop=F]))
  
  if(crve == 'CR0'){
    A.h = Map(function(x) (diag(ncol(x))), x=Xdw.h)
  }else if(crve == 'CR2'){
    A.h = adjustment.matrices.NAAMW(Xdw.h, XTX.inv, delta=-0.5)
  }else if(crve == 'CR3'){
    A.h = adjustment.matrices.NAAMW(Xdw.h, XTX.inv, delta=-1)
  }else{
    stop("p.value.meis() error: wrong crve.type")
  }
  
  model.0 = list()
  model.0$Xdw.h = Xdw.h
  model.0$XTX.inv = XTX.inv
  model.0$A.h = A.h
  model.0$clu.grps = clu.grps
  
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
    
    ind.lo = clu.grps$id.min
    ind.hi = clu.grps$id.max
    nclu = nrow(clu.grps)
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
  # m1 = meis.test.lambdas.1(model.0$Xdw.h, model.0$XTX.inv, model.0$A.h, c2)
  # p0 = p.meis(m1, c2, t1)
  lambdas = bm.test.lambdas(model.0$Xdw.h, model.0$XTX.inv, model.0$A.h, c2, atype='NAAMW')
  p0 = p.meis(model.0$XTX.inv, c2, t1, lambdas=lambdas)
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
    stop(paste(calling.function, 'error: y.name not in data'))
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
  VAR_0 = Map(function(x, eh) (t(a) %*% (t(x) %*% eh) %*% (t(eh) %*% x) %*% a),
              x=Xdw.h, eh=Eh.h, a=model.0$A.h)
  VAR_1 = Reduce('+', VAR_0)
  VAR_2 = model.0$XTX.inv %*% VAR_1 %*% model.0$XTX.inv
  
  #
  # pass inputs to other functions
  #
  c2 = t(t(c1))
  # m1 = meis.test.lambdas.1(model.0$Xdw.h, model.0$XTX.inv, model.0$A.h, c2)
  # q.star = q.meis(m1, c2, alpha=alpha)
  
  lambdas = bm.test.lambdas(model.0$Xdw.h, model.0$XTX.inv, model.0$A.h, c2, atype='NAAMW')
  q.star = q.meis(model.0$XTX.inv, c2, alpha=alpha, lambdas=lambdas)
  
  
  #
  # Construct confidence interval
  #
  lb = (t(c2) %*% B.hat) - (q.star * (t(c2) %*% VAR_2 %*% c2))
  ub = (t(c2) %*% B.hat) + (q.star * (t(c2) %*% VAR_2 %*% c2))
  ci.0 = c(lower.bound=lb, upper.bound=ub)
  
  return(ci.0)
}


#' OLS with a standard battery of tests
#'
#' This function calculates OLS, gets cluster-robust standard errors (CR0),
#' and tests a linear hypothesis for each covariate three ways:
#' assuming the t-stat is normal, assuming the t-stat is t-distributed with
#' degrees of freedom equal to the effective number of clusters, and
#' according to Meiselman (2021).
#'
#' @param data Data frame or matrix containing all outcomes, covariates, and group ids
#' @param y.name Name of dependent variable
#' @param x.names Vector of covariate names
#' @param cluster.id Name of cluster id variable
#' @param weights Name of weight variable, if applicable
#' @param fe.id Name of fixed effect group id variable (defaults to cluster.id)
#' @return numeric, the lower and upper bound of the confidence interval
ols.with.the.works <- function(data, y.name, x.names, cluster.id, weights=NULL, fe.id=cluster.id){
  # 
  # reorder data
  # 
  model.0 = absorb.reorg(data=data, x.names=x.names,
                         cluster.id=cluster.id, fe.id=fe.id, weights=weights, crve='CR0', y.name=y.name)
  
  data.0 = data
  if(!is.null(weights)){
    data.0$reg.weights = data.0[, weights]
  }else{
    data.0$reg.weights = 1
  }
  W = as.matrix(data.0[, 'reg.weights', drop=F])
  Y = as.matrix(data.0[, c(y.name), drop=F])
  X = as.matrix(data.0[, c(x.names), drop=F])
  Xdw.h = model.0$Xdw.h
  Ydw.h = model.0$Ydw.h
  Xdw = do.call(rbind, Xdw.h)
  Ydw = do.call(rbind, Ydw.h)
  
  
  XTY = t(Xdw) %*% Ydw
  Bhat = model.0$XTX.inv %*% XTY
  Ehat = Y - (X %*% Bhat)
  Ew = Ehat * sqrt(W)
  
  #
  # Cluster-robust standard errors
  # 
  
  ind.lo = model.0$clu.grps$id.min
  ind.hi = model.0$clu.grps$id.max
  ngrp = nrow(model.0$clu.grps)
  Ew.h = lapply(seq_len(ngrp), function(x) as.matrix(Ew[ind.lo[x]:ind.hi[x],, drop=F]))
  
  XeeX.h = Map(function(x, e){
    xe = t(x) %*% e
    return(xe %*% t(xe))
  }, x=Xdw.h, e=Ew.h)
  XeeX = Reduce('+', XeeX.h)
  V = model.0$XTX.inv %*% XeeX %*% model.0$XTX.inv
  SE = sqrt(diag(V))
  
  coeff.0 = data.frame(b=c(Bhat), se.cr0=c(SE), t.stat=c(Bhat)/c(SE))
  coeff.0$se.cr1 = coeff.0$se.cr0 * sqrt(ngrp/(ngrp - 1))
  
  #
  # Effective number of clusters
  # for one linear hypothesis on each covariate
  #
  
  c0 = rep(0, length(x.names))
  for(i in 1:length(x.names)){
    c1 = c0
    c1[i] = 1
    c2 = t(t(c1))
    gam.0 = css.test.gammas(Xdw.h, model.0$XTX.inv, c2)
    g.star = (sum(gam.0) ^ 2) / (sum(gam.0 ^ 2))
    coeff.0[i, 'g.star'] = g.star
    
    t0 = coeff.0[i, 't.stat']
    # m1 = meis.test.lambdas.1(model.0$Xdw.h, model.0$XTX.inv, model.0$A.h, c2)
    
    lambdas = bm.test.lambdas(model.0$Xdw.h, model.0$XTX.inv, model.0$A.h, c2, atype='NAAMW')
    df.bm = (sum(lambdas) ^ 2) / (sum(lambdas ^ 2))
    coeff.0[i, 'df.bm'] = df.bm
    
    coeff.0[i, 'p.norm'] = (1 - pnorm(abs(t0))) / 2
    coeff.0[i, 'p.css'] = (1 - pt(t0, df=g.star)) / 2
    coeff.0[i, 'p.meis'] = p.meis(model.0$XTX.inv, c2, t0, lambdas=lambdas)
    
  }
  
  return(coeff.0)
}
