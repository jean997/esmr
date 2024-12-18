
# Calculates E_q[log lik(l, f | Y, omega)]
# The objective function is E_q[log lik(l, f | Y, omega)] - KL(q_l || g_l)
## E[ sum (Y, - y_j)^T Omega (Y_j - y_j)]
## This version treats all betas as independent
#'@export
calc_ell2 <- function(Y, abar, a2bar, fgbar, omega, omega_logdet, s_equal){
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
    ell <- -0.5*(part_a + part_b - omega_logdet)
  }else{
    part_a <- map_dbl(seq(n), function(i){
      r <- R[i,,drop = FALSE]
      tcrossprod(r, tcrossprod(r, omega[[i]]))
    }) %>% sum()

    diagA <- map(seq(n), function(i){
      colSums(crossprod(omega[[i]],fgbar) * fgbar)
    }) %>% unlist() %>% matrix(ncol = k, byrow = T)

    part_b <- sum(varabar*diagA)

    ell <- -0.5*(part_a + part_b  - omega_logdet)
  }

  return(ell)
}


#'@export
logLik.esmr <- function(x, ...) {
  log_py(x, ...)
}








