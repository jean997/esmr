

#'@export
update_lj_new <- function(dat, j,
                          g_init = NULL, fix_g = FALSE,
                          return_post = FALSE, return_x_s = FALSE){

  R_j <- dat$Y - (dat$l$abar[,-j,drop=FALSE] %*% t(dat$f$fgbar[,-j,drop=FALSE]))
  fgbar_j <- with(dat$f, fgbar[,j])
  HxH <- outer(fgbar_j, fgbar_j)
  g_j <- dat$G[,j]
  gA_j <- dat$G[,-j, drop = FALSE] %*% t(dat$l$abar[,-j, drop = FALSE])

  if(dat$s_equal){
    a1 <- (HxH * dat$omega) %>% sum()
    a2 <- crossprod(crossprod(dat$f$fVo, g_j), g_j)
    A <- rep(a1 + a2, n)
    b1 <- R_j %*% dat$omega %*% fgbar_j
    b2 <- g_j %*% dat$f$fVo %*% gA_j
    B <- b1-b2
  }else{
    A <- map(seq(dat$n), function(i){
      a1 <- sum(HxH * dat$omega[[i]])
      a2 <- crossprod(crossprod(dat$f$fVo[[i]], g_j), g_j) #t(g_j) %*% dat$beta$VVo[[i]] %*% g_j
      a1 + a2
    }) %>% unlist()
    B <- map(seq(dat$n), function(i){
      b1 <- R_j[i,, drop = FALSE] %*% dat$omega[[i]] %*% fgbar_j
      b2 <- g_j %*% dat$f$fVo[[i]] %*% gA_j[,i,drop = FALSE]
      b1 -b2
    }) %>% unlist()
  }

  x <- B/A
  s <- 1/sqrt(A)

  if(return_x_s){
    return(list(x = x, s = s))
  }

  ixnmiss <- which(A > 0)
  if(length(ixnmiss)  != dat$n){
    x <- x[ixnmiss]
    s <- s[ixnmiss]
  }


  ebnm_res <- dat$ebnm_fn( x= as.numeric(x), s = s, g_init = g_init, fix_g= fix_g, output = ebnm::ebnm_output_all())
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
