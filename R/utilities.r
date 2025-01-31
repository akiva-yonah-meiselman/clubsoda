

# library('pracma')

#' Matrix Exponent
#'
#' This function raises a matrix to an exponent using
#' the matrix's eigendecomposition.
#'
#' @param x A square matrix
#' @param p Linear hypothesis vector
#' @param symmetric Indicates whether x is a symmetric matrix (if true, improves performance slightly)
#' @param tol Error tolerance for eigenvalues close to zero assumed to be zero
#' @return A vector of eigenvalues
matrix_power <- function(x, p, symmetric = TRUE, tol = -12) {
  eig <- eigen(x, symmetric = symmetric)
  v1 <- eig$values / max(abs(eig$values))
  v2 <- ifelse(abs(v1) > 10^tol, v1^p, 0)
  v3 <- v2 * (max(abs(eig$values)) ^ p)
  with(eig, vectors %*% (v3 * t(vectors)))
}

# 
# CRVE adjustment matrices
# 

#' Adjustment Matrices
#'
#' This function gets adjustment matrices for the
#' cluster-robust variance estimators CR2 and CR3
#' from Bell and McCaffrey (2002). The adjustment
#' matrices are given by Niccodemi et al (2020)
#' to be calculated much faster than those in BM.
#' @param xg xg List of cluster-specific submatrices of design matrix
#' @param xtx Inverse of design matrix squared
#' @param delta Power to raise (I_g - H_gg), by default delta=-0.5 for CR2
#' @return A vector of eigenvalues, first not yet scaled by q
adjustment.matrices.NAAMW <- function(xg, xtx, delta=-0.5, xgtxg=NULL){
  if(is.null(xgtxg)){
    xgtxg = Map(function(x) (t(x) %*% x), x=xg)
  }
  
  xtx.inv.half = clubsoda:::matrix_power(xtx, p=0.5)
  xtx.half = solve(xtx.inv.half)
  ik = diag(ncol(xtx))
  ag = Map(function(x){
    meat.0 = (ik - ( xtx.inv.half %*% x %*% xtx.inv.half ))
    meat.1 = clubsoda:::matrix_power(meat.0, p=delta)
    return(xtx.inv.half %*% meat.1 %*% xtx.half)
  }, x=xgtxg)
  return(ag)
}

#' Adjustment Matrices
#'
#' This function gets adjustment matrices for the
#' cluster-robust variance estimators CR2 and CR3
#' from Bell and McCaffrey (2002). The is a
#' straightforward calculation of A_g.
#' 
#' @param xg List of cluster-specific submatrices of design matrix
#' @param xtx Inverse of design matrix squared
#' @param delta Power to raise (I_g - H_gg), by default delta=0.5 for CR2
#' @return A vector of eigenvalues
adjustment.matrices.BM <- function(xg, xtx, delta=-0.5) {
  a.g = Map(function(x){
    I_g = diag(nrow(x))
    H_gg = x %*% xtx %*% t(x)
    a = matrix_power(I_g - H_gg, p=delta)
    return(a)
  }, x=xg)
  return(a.g)
}





