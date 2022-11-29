update_l_sequential <- function(dat, jj){
  l_update <- list()
  lbar <- dat$l$lbar
  l2bar <- dat$l$l2bar

  wpost <- dat$l$wpost
  mupost <- dat$l$mupost
  s2post <- dat$l$s2post


  lfsr <- dat$l$lfsr
  g_hat <- dat$l$g_hat
  if(!missing(jj)){
    coords <- jj
  }else{
    coords <- seq(dat$p)
  }

  for(j in coords){
    R_j <- dat$Y - (lbar[,-j,drop=FALSE] %*% t(dat$f$fbar[,-j,drop=FALSE]))
    lu <- update_l_k(R_j, dat$f$fbar[,j], dat$f$f2bar[,j], dat$omega, dat$ebnm_fn)

    lbar[lu$posterior$index,j] <- lu$posterior$mean
    l2bar[lu$posterior$index,j] <- lu$posterior$second_moment
    wpost[lu$posterior$index,j] <- lu$posterior$wpost
    mupost[lu$posterior$index,j] <- lu$posterior$mu
    s2post[lu$posterior$index,j] <- lu$posterior$s2
    lfsr[lu$posterior$index,j] <- lu$posterior$lfsr
    g_hat[[j]] <- lu$fitted_g
    l_update[[j]] <- lu
  }
  kl <- map(l_update, "KL") %>% unlist() %>% sum()
  dat$l <- list(lbar =lbar, l2bar = l2bar, lfsr = lfsr,
                wpost = wpost, mupost = mupost, s2post = s2post,
                kl = kl, g_hat = g_hat)
  return(dat)
}


update_beta_sequential <- function(dat, jj){
  p <- dat$p

  beta_j <- dat$beta$beta_j
  beta_k <- dat$beta$beta_k

  fbar <- dat$f$fbar
  f2bar <- dat$f$f2bar
  beta_m <- dat$beta$beta_m
  beta_s <- dat$beta$beta_s

  if(missing(jj)){
    coords <- seq(length(beta_j))
  }else{
    coords <- jj
  }


  for(i in coords){
      k <- beta_k[i]
      j <- beta_j[i]

      R_k <- dat$Y - (dat$l$lbar[,-k,drop=FALSE]%*%t(fbar[,-k,drop=FALSE]))
      b <- update_beta_k(R_k = R_k, j=j, k=k,
                    lbar=dat$l$lbar, l2bar=dat$l$l2bar,
                    omega = dat$omega, fbar = fbar,
                    sigma_beta = dat$sigma_beta)
      beta_m[i] <- b$m
      beta_s[i] <- b$s
      f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
      fbar <- f$fbar
      f2bar <- f$f2bar
  }
  dat$beta$beta_m <- beta_m
  dat$beta$beta_s <- beta_s
  dat$f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
  return(dat)
}

update_l_sequential_future <- function(dat, jj){
  l_update <- list()
  lbar_o <- dat$l$lbar_o
  l2bar_o <- dat$l$l2bar_o

  wpost <- dat$l$wpost
  mupost <- dat$l$mupost
  s2post <- dat$l$s2post


  lfsr <- dat$l$lfsr
  g_hat <- dat$l$g_hat
  if(!missing(jj)){
    coords <- jj
  }else{
    coords <- seq(dat$p)
  }

  for(j in coords){
    R_j <- dat$Y - (lbar_o[,-j,drop=FALSE] %*% t(dat$f$fbar[,-j,drop=FALSE]))
    lu <- update_l_k(R_j, dat$f$fbar[,j], dat$f$f2bar[,j], dat$omega, dat$ebnm_fn)

    lbar_o[lu$posterior$index,j] <- lu$posterior$mean
    l2bar_o[lu$posterior$index,j] <- lu$posterior$second_moment
    wpost[lu$posterior$index,j] <- lu$posterior$wpost
    mupost[lu$posterior$index,j] <- lu$posterior$mu
    s2post[lu$posterior$index,j] <- lu$posterior$s2
    lfsr[lu$posterior$index,j] <- lu$posterior$lfsr
    g_hat[[j]] <- lu$fitted_g
    l_update[[j]] <- lu
  }
  kl <- map(l_update, "KL") %>% unlist() %>% sum()


  lbar <- lbar_o %*% t(dat$G)
  Vl <- l2bar_o - (lbar_o^2)
  l2bar <- (lbar^2) + (Vl %*% t(dat$G)^2)

  dat$l <- list(lbar =lbar, l2bar = l2bar,
                lbar_o = lbar_o, l2bar_o= l2bar_o,
                lfsr = lfsr,
                wpost = wpost, mupost = mupost, s2post = s2post,
                kl = kl, g_hat = g_hat)
  return(dat)
}




update_beta_sequential_future <- function(dat, jj){
  p <- dat$p

  beta_j <- dat$beta$beta_j
  beta_k <- dat$beta$beta_k

  fbar <- dat$f$fbar_o
  f2bar <- dat$f$f2bar_o
  beta_m <- dat$beta$beta_m
  beta_s <- dat$beta$beta_s

  if(missing(jj)){
    coords <- seq(length(beta_j))
  }else{
    coords <- jj
  }


  for(i in coords){
    k <- beta_k[i]
    j <- beta_j[i]

    R_k <- dat$Y - dat$l$lbar[,-k,drop=FALSE]%*%t(fbar[,-k,drop=FALSE])
    b <- update_beta_k(R_k = R_k, j=j, k=k,
                       lbar=dat$l$lbar, l2bar=dat$l$l2bar,
                       omega = dat$omega, fbar = fbar,
                       sigma_beta = dat$sigma_beta)
    beta_m[i] <- b$m
    beta_s[i] <- b$s
    f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
    fbar <- f$fbar_o
    f2bar <- f$f2bar_o
  }
  dat$beta$beta_m <- beta_m
  dat$beta$beta_s <- beta_s
  dat$f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
  return(dat)
}

