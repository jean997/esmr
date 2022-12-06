
#'@export
update_l_k <- function(R_k, fbar_k, f2bar_k, omega, ebnm_fn){
  n <- nrow(R_k)
  p <- ncol(R_k)
  stopifnot(length(fbar_k) == p)
  stopifnot(length(f2bar_k) == p)
  s_equal <- check_equal_omega(omega)


  fbar_k <- matrix(fbar_k, ncol = 1)
  Fbar <- fbar_k %*% t(fbar_k)
  #diag(Fbar) <- f2bar_k
  # Fbar is 2 by 2
  if(s_equal){
    A <- (Fbar * omega) %>% sum()
    A <- rep(A, n)
    B <- R_k %*% omega %*% fbar_k
  }else{
    A <- map(omega, function(o){
      (Fbar * o) %>% sum()
    }) %>% unlist()
    B <- map(seq(n), function(i){
      R_i <- R_k[i,, drop = FALSE]
      R_i %*% omega[[i]] %*% fbar_k
    }) %>% unlist()
  }

  x <- B/A
  s <- 1/sqrt(A)

  ixnmiss <- which(A > 0)
  if(length(ixnmiss)  != n){
    x <- x[ixnmiss]
    s <- s[ixnmiss]
  }


  ebnm_res <- ebnm_fn( x= as.numeric(x), s = s, g_init = NULL, fix_g= FALSE, output = output_all())
  ebnm_res$KL <-  (ebnm_res$log_likelihood
                   - flashier:::normal.means.loglik(x,s,
                                                    ebnm_res$posterior$mean,
                                                    ebnm_res$posterior$second_moment))
  ebnm_res$posterior$index <- ixnmiss
  # This is only for point normal
  # a <- 1/ebnm_res$fitted_g$sd[2]^2
  # w <- ebnm_res$fitted_g$pi[2]
  # ebnm_res$posterior$wpost <- ebnm:::wpost_normal(x, s, w, a, 0)
  # ebnm_res$posterior$mu <- ebnm:::pmean_cond_normal(x, s, a, 0)
  # ebnm_res$posterior$s2 <- ebnm:::pvar_cond_normal(s, a)
  return(ebnm_res)


}

# ebnm_res$KL  is computed as log p(x | g) - E_{p(theta | g)}[p(x | theta)] = E_{p(theta | g)}(log p(theta) - log(p(theta | x)))
# = -KL(p(theta | x) || p(theta))
