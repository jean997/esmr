update_l_sequential_2parts <- function(dat, jj){
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

  ix1 <- dat$ix1
  ix0 <- dat$ix0

  for(j in coords){
    R_j <- matrix(nrow = dat$n, ncol = dat$p)
    R_j[ix0,] <- dat$Y[ix0,] - (abar[ix0,-j,drop=FALSE] %*% t(dat$f$fgbar[,-j,drop=FALSE]))
    R_j[ix1,] <- dat$Y[ix1,] - (abar[ix1,-j,drop=FALSE] %*% t(dat$f0$fgbar[,-j,drop=FALSE]))
    lu <- update_l_k_2parts(R_j, dat$f$fgbar[,j], dat$f$fg2bar[,j],
                            dat$f0$fgbar[,j], dat$f0$fg2bar[,j],
                            dat$omega, dat$ebnm_fn, ix0)

    abar[lu$posterior$index,j] <- lu$posterior$mean
    a2bar[lu$posterior$index,j] <- lu$posterior$second_moment

    lfsr[lu$posterior$index,j] <- lu$posterior$lfsr
    g_hat[[j]] <- lu$fitted_g
    l_update[[j]] <- lu
  }
  kl <- map(l_update, "KL") %>% unlist() %>% sum()

  lbar <- abar %*% t(dat$G)
  Va <- a2bar - (abar^2)
  l2bar <- (lbar^2) + (Va %*% t(dat$G)^2)

  dat$l <- list(lbar =lbar, l2bar = l2bar,
                abar = abar, a2bar= a2bar,
                lfsr = lfsr,
                kl = kl, g_hat = g_hat)
  return(dat)
}






#'@export
update_l_k_2parts <- function(R_k, fgbar_k, fg2bar_k,
                              fgbar0_k, fg2bar0_k,
                              omega, ebnm_fn,
                              ix0){
  n <- nrow(R_k)
  p <- ncol(R_k)
  stopifnot(length(fgbar_k) == p)
  stopifnot(length(fg2bar_k) == p)
  s_equal <- check_equal_omega(omega)


  fgbar_k <- matrix(fgbar_k, ncol = 1)
  fgbar <- fgbar_k %*% t(fgbar_k)
  diag(fgbar) <- fg2bar_k

  fgbar0_k <- matrix(fgbar0_k, ncol = 1)
  fgbar0 <- fgbar0_k %*% t(fgbar0_k)
  diag(fgbar0) <- fg2bar0_k

  # fgbar is 2 by 2
  if(s_equal){
    A1 <- (fgbar * omega) %>% sum()
    A0 <- (fgbar0 * omega) %>% sum()
    A <- rep(A1, n)
    A[ix0] <- A0
    B1 <- R_k %*% omega %*% fgbar_k
    B0 <- R_k %*% omega %*% fgbar0_k
    B <- B1
    B[ix0] <- B0[ix0]
  }else{
    A1 <- map(omega, function(o){
      (fgbar * o) %>% sum()
    }) %>% unlist()
    A0 <- map(omega, function(o){
      (fgbar0 * o) %>% sum()
    }) %>% unlist()
    A <- A1
    A[ix0] <- A0[ix0]
    B1 <- map(seq(n), function(i){
      R_i <- R_k[i,, drop = FALSE]
      R_i %*% omega[[i]] %*% fgbar_k
    }) %>% unlist()
    B0 <- map(seq(n), function(i){
      R_i <- R_k[i,, drop = FALSE]
      R_i %*% omega[[i]] %*% fgbar0_k
    }) %>% unlist()
    B <- B1
    B[ix0] <- B0[ix0]
  }

  x <- B/A
  s <- 1/sqrt(A)

  ixnmiss <- which(A > 0)
  if(length(ixnmiss)  != n){
    x <- x[ixnmiss]
    s <- s[ixnmiss]
  }


  ebnm_res <- ebnm_fn( x= as.numeric(x), s = s, g_init = NULL, fix_g= FALSE, output = ebnm::ebnm_output_all())
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
  # ebnm_res$posterior$post_mode <- round(ebnm_res$posterior$wpost)*ebnm_res$posterior$mu
  return(ebnm_res)


}

# ebnm_res$KL  is computed as log p(x | g) - E_{p(theta | g)}[p(x | theta)] = E_{p(theta | g)}(log p(theta) - log(p(theta | x)))
# = -KL(p(theta | x) || p(theta))
