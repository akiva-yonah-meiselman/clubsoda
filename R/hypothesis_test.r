


# library('CompQuadForm')


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
#' @export
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
#' @export
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
#' @export
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
#' @export
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





