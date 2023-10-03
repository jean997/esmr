
# Calculates E_q[log lik(l, f | Y, omega)]
# The objective function is E_q[log lik(l, f | Y, omega)] - KL(q_l || g_l)
## E[ sum (Y, - y_j)^T Omega (Y_j - y_j)]
#'@export
calc_ell2 <- function(Y, abar, a2bar, fgbar, omega){
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(fgbar)
  check_matrix(abar, "abar", n, k)
  check_matrix(a2bar, "a2bar", n, k)
  check_matrix(fgbar, "fgbar", p, k)

  s_equal <- check_equal_omega(omega)

  ybar <- abar %*% t(fgbar)
  R <- Y - ybar

  varabar <- a2bar - (abar^2)

  if(s_equal){
    part_a <- sum(tcrossprod(R, omega) * R) #quad.tdiag(omega, R) %>% sum()
    diagA <- colSums(crossprod(omega,fgbar) * fgbar) #quad.diag(omega, fgbar) # diag( t(F) %*% omega %*% F)
    part_b <- t(t(varabar)*diagA) %>% sum()
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

    part_b <- sum(varabar*diagA)
    #part_b <- t(t(varlbar)*diagA) %>% sum()
    ell <- -0.5*(part_a + part_b)
  }
  return(ell)
}


## calculate E_q[log g(a_i) - log q(a_i)]
## E_q[log g(a_i)] = E_q[ log( \sum_K pi_k P(a_i; m_k, s_k))]
calc_kl_k <- function(x, s, w, a, mu, post_mean, post_second){
  #n <- length(x)
  part1 <- ebnm:::loglik_point_normal(x, s, w, a, mu)
  part2 <- flashier:::normal.means.loglik(x,s,
                                          post_mean,
                                          post_second)
  # wpost <- ebnm_res$posterior$wpost
  # apost <- 1/ebnm_res$posterior$s2
  # mupost <- ebnm_res$posterior$mu
  # part3 <- sapply(1:n, function(i){
  #   ebnm:::loglik_point_normal(x[i], s[i], wpost[i], apost[i], mupost[i] )
  # }) |> sum()

  return(part1 - part2)
}
