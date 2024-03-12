update_l_sequential <- function(dat, jj){
  l_update <- list()
  abar <- dat$l$abar
  a2bar <- dat$l$a2bar

  lfsr <- dat$l$lfsr
  g_hat <- dat$l$g_hat
  if(!missing(jj)){
    coords <- jj
  }else{
    coords <- seq(dat$k)
  }

  for(j in coords){
    R_j <- dat$Y - (abar[,-j,drop=FALSE] %*% t(dat$f$fgbar[,-j,drop=FALSE]))
    lu <- update_l_k(R_j, dat$f$fgbar[,j], dat$f$fg2bar[,j], dat$omega, dat$ebnm_fn)

    abar[lu$posterior$index,j] <- lu$posterior$mean
    a2bar[lu$posterior$index,j] <- lu$posterior$second_moment

    lfsr[lu$posterior$index,j] <- lu$posterior$lfsr
    g_hat[[j]] <- lu$fitted_g
    l_update[[j]] <- lu
  }
  kl <- map(l_update, "KL") %>% unlist() %>% sum()

  lbar <- abar %*% t(dat$G)
  all_zero_cols_lbar <- apply(zapsmall(lbar, digits = 10), 2, function(x){
    all(x == 0)
  })
  if (any(all_zero_cols_lbar)) {
    stop('lbar has a column of all zeros for column(s): ', which(all_zero_cols_lbar))
  }
  Va <- a2bar - (abar^2)
  l2bar <- (lbar^2) + (Va %*% t(dat$G)^2)

  dat$l <- list(lbar =lbar, l2bar = l2bar,
                abar = abar, a2bar= a2bar,
                lfsr = lfsr,
                kl = kl, g_hat = g_hat)
  return(dat)
}






#'@export
update_l_k <- function(R_k, fgbar_k, fg2bar_k, omega, ebnm_fn,
                       g_init = NULL, fix_g = FALSE,
                       return_post = FALSE, return_x_s = FALSE){
  n <- nrow(R_k)
  p <- ncol(R_k)
  stopifnot(length(fgbar_k) == p)
  stopifnot(length(fg2bar_k) == p)
  s_equal <- check_equal_omega(omega)


  fgbar_k <- matrix(fgbar_k, ncol = 1)
  fgbar <- fgbar_k %*% t(fgbar_k)
  diag(fgbar) <- fg2bar_k
  # fgbar is 2 by 2
  if(s_equal){
    A <- (fgbar * omega) %>% sum()
    A <- rep(A, n)
    B <- R_k %*% omega %*% fgbar_k
  }else{
    A <- map(omega, function(o){
      (fgbar * o) %>% sum()
    }) %>% unlist()
    B <- map(seq(n), function(i){
      R_i <- R_k[i,, drop = FALSE]
      R_i %*% omega[[i]] %*% fgbar_k
    }) %>% unlist()
  }

  x <- B/A
  s <- 1/sqrt(A)

  if(return_x_s){
    return(list(x = x, s = s))
  }

  ixnmiss <- which(A > 0)
  if(length(ixnmiss)  != n){
    x <- x[ixnmiss]
    s <- s[ixnmiss]
  }


  ebnm_res <- ebnm_fn( x= as.numeric(x), s = s, g_init = g_init, fix_g= fix_g, output = ebnm::ebnm_output_all())
  ebnm_res$KL <-  (ebnm_res$log_likelihood
                   - flashier:::normal.means.loglik(x,s,
                                                    ebnm_res$posterior$mean,
                                                    ebnm_res$posterior$second_moment))
  ebnm_res$posterior$index <- ixnmiss
  # This is only for point normal
  if(return_post){
    a <- 1/ebnm_res$fitted_g$sd[2]^2
    w <- ebnm_res$fitted_g$pi[2]
    ebnm_res$posterior$wpost <- ebnm:::wpost_normal(x, s, w, a, 0)
    ebnm_res$posterior$mu <- ebnm:::pmean_cond_normal(x, s, a, 0)
    ebnm_res$posterior$s2 <- ebnm:::pvar_cond_normal(s, a)
    #ebnm_res$posterior$post_mode <- round(ebnm_res$posterior$wpost)*ebnm_res$posterior$mu
  }
  return(ebnm_res)


}

# ebnm_res$KL  is computed as log p(x | g) - E_{p(theta | g)}[p(x | theta)] = E_{p(theta | g)}(log p(theta) - log(p(theta | x)))
# = -KL(p(theta | x) || p(theta))
