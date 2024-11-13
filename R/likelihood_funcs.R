
# Calculates E_q[log lik(l, f | Y, omega)]
# The objective function is E_q[log lik(l, f | Y, omega)] - KL(q_l || g_l)
## E[ sum (Y, - y_j)^T Omega (Y_j - y_j)]
## This version treats all betas as independent
#'@export
calc_ell2 <- function(Y, abar, a2bar, fgbar, omega, s_equal, R_is_id = FALSE){
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(fgbar)

  #fgbar <- fbar %*% G

  #s_equal <- check_equal_omega(omega)

  ybar <- abar %*% t(fgbar)
  R <- Y - ybar
  varabar <- a2bar - (abar^2)

  if(s_equal){
    part_a <- sum(tcrossprod(R, omega) * R) #quad.tdiag(omega, R) %>% sum()
    diagA <- colSums(crossprod(omega,fgbar) * fgbar)
    part_b <- t(t(varabar)*diagA) %>% sum()
    # TODO: For now if omega is diagonal then we can just sum the log diagonal
    # Might need change the determinant part for larger graphs too
    if (R_is_id) {
      log_det <- n * sum(log(diag(omega)))
    } else {
      log_det <- log(det(omega))*n
    }
    ell <- -0.5*(part_a + part_b -log_det)
  }else{
    part_a <- map_dbl(seq(n), function(i){
      r <- R[i,,drop = FALSE]
      tcrossprod(r, tcrossprod(r, omega[[i]]))
    }) %>% sum()

    diagA <- map(seq(n), function(i){
      colSums(crossprod(omega[[i]],fgbar) * fgbar)
    }) %>% unlist() %>% matrix(ncol = k, byrow = T)

    part_b <- sum(varabar*diagA)

    if (R_is_id) {
      log_det <- sapply(omega, function(o){sum(log(diag(o)))*n}) %>% sum()
    } else {
      log_det <- sapply(omega, function(o){log(det(o))}) %>% sum()
    }
    ell <- -0.5*(part_a + part_b  - log_det)
  }

  return(ell)
}


#'@export
logLik.esmr <- function(x) {
  log_py(x)
}








