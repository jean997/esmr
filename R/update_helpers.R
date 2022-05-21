update_l_sequential <- function(dat, ebnm_fn){
  l_update <- list()
  lbar <- dat$l$lbar
  l2bar <- dat$l$l2bar
  for(j in seq(p)){
    R_j <- dat$Y - (lbar[,-j,drop=FALSE] %*% t(dat$f$fbar[,-j,drop=FALSE]))
    lu <- update_l_k(R_j, dat$f$fbar[,j], dat$f$f2bar[,j], dat$omega, ebnm_fn)

    lbar[lu$posterior$index,j] <- lu$posterior$mean
    l2bar[lu$posterior$index,j] <- lu$posterior$second_moment
    l_update[[j]] <- lu
  }
  kl <- map(l_update, "KL") %>% unlist() %>% sum()
  dat$l <- list(lbar =lbar, l2bar = l2bar, kl = kl)
  return(dat)
}


update_beta_sequential <- function(dat, sigma_beta){
  p <- dat$p

  beta_j <- dat$beta$beta_j
  beta_k <- dat$beta$beta_k
  ix <- dat$ix_beta
  fbar <- dat$f$fbar
  f2bar <- dat$f$f2bar
  beta_m <- dat$beta$beta_m
  beta_s <- dat$beta$beta_s

  for(i in seq(length(beta_j))){
      k <- beta_k[i]
      j <- beta_j[i]
      R_k <- dat$Y[ix,] - (dat$l$lbar[ix,-k,drop=FALSE]%*%t(fbar[,-k,drop=FALSE]))
      b <- update_beta_k(R_k = R_k, j=j, k=k,
                    lbar=dat$l$lbar[ix,], l2bar=dat$l$l2bar[ix,],
                    omega = dat$omega, fbar = fbar,
                    sigma_beta = sigma_beta)
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

  check <- 1
  obj <- c()
  obj_old <- -Inf
  i <- 1
  while(i < max_iter & check > tol){

    dat <- update_l_sequential(dat, dat$ebnm_fn)


    # beta update
    if(!fix_beta & !beta_joint){
      dat <- update_beta_sequential(dat, dat$sigma_beta)

    }else if(!fix_beta & beta_joint){
      beta_upd <- with(dat,
                       update_beta_joint(Y[ix_beta,], l$lbar[ix_beta,], l$l2bar[ix_beta,], omega))
      dat$beta$beta_m <- beta_upd$m
      dat$beta$beta_s <- sqrt(diag(beta_upd$S))
      dat$beta$beta_var <- beta_upd$S
      dat$f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
    }



    obj <- c(obj, dat$l$kl)
    obj_new <- obj[length(obj)]
    check <- obj_new - obj_old
    obj_old <- obj_new

    check <- abs(check)
    cat(i, ": ", obj_new, " ", dat$beta$beta_m, "\n")
    i <- i + 1
  }
  dat$obj <- obj
  return(dat)
}
