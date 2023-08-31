
# Calculates E_q[log lik(l, f | Y, omega)]
# The objective function is E_q[log lik(l, f | Y, omega)] - KL(q_l || g_l)
#'@export
calc_ell2 <- function(Y, lbar, l2bar, fgbar, omega){
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(fgbar)
  check_matrix(lbar, "lbar", n, k)
  check_matrix(l2bar, "l2bar", n, k)
  check_matrix(fgbar, "fgbar", p, k)

  s_equal <- check_equal_omega(omega)

  ybar <- lbar %*% t(fgbar)
  R <- Y - ybar

  varlbar <- l2bar - (lbar^2)

  if(s_equal){
    part_a <- sum(tcrossprod(R, omega) * R) #quad.tdiag(omega, R) %>% sum()
    diagA <- colSums(crossprod(omega,fgbar) * fgbar) #quad.diag(omega, fgbar) # diag( t(F) %*% omega %*% F)
    part_b <- t(t(varlbar)*diagA) %>% sum()
    ell <- -0.5*(part_a + part_b)
  }else{
    part_a <- map_dbl(seq(n), function(i){
      #quad.tform(omega[[i]], R[i,,drop = FALSE])
      r <- R[i,,drop = FALSE]
      tcrossprod(r, tcrossprod(r, omega[[i]]))
    }) %>% sum()

    diagA <- map(seq(n), function(i){
      colSums(crossprod(omega[[i]],fgbar) * fgbar)
    }) %>% unlist() %>% matrix(ncol = k, byrow = T)

    part_b <- sum(varlbar*diagA)
    #part_b <- t(t(varlbar)*diagA) %>% sum()
    ell <- -0.5*(part_a + part_b)
  }
  return(ell)
}
