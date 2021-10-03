

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
#' from Bell and McCaffrey (2002).
#' 
#' @param xg List of cluster-specific submatrices of design matrix
#' @param xtx Inverse of design matrix squared
#' @param delta Power to raise (I_g - H_gg), by default delta=0.5 for CR2
#' @return A vector of eigenvalues
adjustment.matrices <- function(xg, xtx, delta=0.5) {
  
  xtx.2 = matrix_power(xtx, p=0.5)
  
  a.g = Map(function(x){
    if(ncol(x) < nrow(x)){
      # eigendecompose the inside-out hat submatrix
      h.0 = eigen(xtx.2 %*% (t(x) %*% x) %*% xtx.2)
      h.00 = h.0$values[zapsmall(h.0$values) != 0]
      h.01 = as.matrix(h.0$vectors[, zapsmall(h.0$values) != 0])
      # orthonormal basis for the inside-out hat submatrix
      h.02 = diag(x=h.00 ^ (-0.5), nrow=length(h.00))
      h.1 = x %*% (xtx.2 %*% (h.01) %*% h.02)
      # orthonormal basis for the hat submatrix, H_gg
      h.2 = cbind(h.1, pracma::nullspace(t(h.1)))
      # eigenvalues of (I_g - H_gg)^delta
      h.3 = 1 - c(h.00, rep(0, nrow(x) - length(h.00)))
      h.4 = ifelse(zapsmall(h.3) == 0, 0, h.3 ^ -delta)
      # basis %*% values %*% basis, but faster
      h.5 = sapply(1:ncol(h.2),function(k) h.2[,k] * h.4[k])
      h.6 = h.2 %*% t(h.5)
      return(h.6)
    }else{
      I_g = diag(nrow(x))
      H_gg = x %*% xtx %*% t(x)
      a = matrix_power(I_g - H_gg, p=-delta)
      return(a)
    }
  }, x=xg)
  return(a.g)
}


#' Adjustment Matrices
#'
#' This function gets adjustment matrices for the
#' cluster-robust variance estimators CR2 and CR3
#' from Bell and McCaffrey (2002). The is a more
#' straightforward calculation of adjustment.matrices(),
#' for testing to make sure that adjustment.matrices()
#' is correct.
#' 
#' @param xg List of cluster-specific submatrices of design matrix
#' @param xtx Inverse of design matrix squared
#' @param delta Power to raise (I_g - H_gg), by default delta=0.5 for CR2
#' @return A vector of eigenvalues
adjustment.matrices.slow <- function(xg, xtx, delta=0.5) {
  a.g = Map(function(x){
    I_g = diag(nrow(x))
    H_gg = x %*% xtx %*% t(x)
    a = matrix_power(I_g - H_gg, p=-delta)
    return(a)
  }, x=xg)
  return(a.g)
}





