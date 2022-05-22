update_l_sequential <- function(dat, jj){
  l_update <- list()
  lbar <- dat$l$lbar
  l2bar <- dat$l$l2bar
  lfsr <- dat$l$lfsr
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
    lfsr[lu$posterior$index,j] <- lu$posterior$lfsr
    l_update[[j]] <- lu
  }
  kl <- map(l_update, "KL") %>% unlist() %>% sum()
  dat$l <- list(lbar =lbar, l2bar = l2bar, lfsr = lfsr, kl = kl)
  return(dat)
}


update_beta_sequential <- function(dat, jj){
  p <- dat$p

  beta_j <- dat$beta$beta_j
  beta_k <- dat$beta$beta_k
  #ix <- dat$ix_beta
  fbar <- dat$f$fbar
  f2bar <- dat$f$f2bar
  beta_m <- dat$beta$beta_m
  beta_s <- dat$beta$beta_s

  if(missing(jj)){
    coords <- seq(length(beta_j))
  }else{
    coords <- jj
  }

  ix <- dat$ix

  for(i in coords){
      k <- beta_k[i]
      j <- beta_j[i]
      #ix <- which(dat$lfsr1[,k] < dat$lfsr_thresh)
      #ix <- which(dat$pval[,k] < dat$pval_thresh)

      R_k <- dat$Y[ix,] - (dat$l$lbar[ix,-k,drop=FALSE]%*%t(fbar[,-k,drop=FALSE]))
      b <- update_beta_k(R_k = R_k, j=j, k=k,
                    lbar=dat$l$lbar[ix,], l2bar=dat$l$l2bar[ix,],
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

ebmr_solve <- function(dat, max_iter, tol){

  fix_beta <- dat$fix_beta
  beta_joint <- dat$beta_joint
  est_tau <- dat$est_tau

  check <- 1
  obj <- c()
  obj_old <- -Inf
  i <- 1

  if(dat$est_tau){
    dat$omega_given <- dat$omega
  }
  while(i < max_iter & check > tol){

    dat <- update_l_sequential(dat)
    if(i == 1){
      dat$lfsr1 <- dat$l$lfsr
      dat$lfsr1[,1] <- 0
    }
    #dat$l$lbar[dat$lfsr1 > dat$lfsr_thresh] <- 0
    #dat$l$l2bar[dat$lfsr1 > dat$lfsr_thresh] <- 0
    #dat$l$lbar[dat$pval > dat$pval_thresh] <- 0
    #dat$l$l2bar[dat$pval > dat$pval_thresh] <- 0
    # beta update
    if(!fix_beta & !beta_joint){
      dat <- update_beta_sequential(dat)

    }else if(!fix_beta & beta_joint){
      beta_upd <- with(dat,
                       update_beta_joint(Y[ix,], l$lbar[ix,], l$l2bar[ix,], omega))
      dat$beta$beta_m <- beta_upd$m
      dat$beta$beta_s <- sqrt(diag(beta_upd$S))
      dat$beta$beta_var <- beta_upd$S
      dat$f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
    }
    if(est_tau){
      dat <- update_tau(dat)
    }
    obj <- c(obj, dat$l$kl)
    obj_new <- obj[length(obj)]
    check <- obj_new - obj_old
    obj_old <- obj_new

    check <- abs(check)
    cat(i, ": ", obj_new, " ", dat$beta$beta_m, "\n")
    i <- i + 1
  }
  dat <- update_l_sequential(dat)
  dat$obj <- obj
  return(dat)
}
