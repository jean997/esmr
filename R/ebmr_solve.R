ebmr_solve <- function(dat, max_iter, tol){

  fix_beta <- dat$fix_beta
  beta_joint <- dat$beta_joint
  est_tau <- dat$est_tau

  check <- 1
  obj <- c()
  obj_old <- -Inf
  i <- 1
  if(dat$ll){
    myll <- c()
  }
  if(dat$est_tau){
    dat$omega_given <- dat$omega
  }
  while(i < max_iter & check > tol){

    dat <- update_l_sequential(dat)
    if(dat$ll){
      myll <- c(myll, with(dat, calc_ell2(Y, l$lbar, l$l2bar, f$fbar, omega)))
    }
    if(i == 1){
      dat$lfsr1 <- dat$l$lfsr
      dat$lfsr1[,1] <- 0
    }
    # beta update
    if(!fix_beta & !beta_joint){
      dat <- update_beta_sequential(dat)

    }else if(!fix_beta & beta_joint){
      beta_upd <- with(dat,
                       update_beta_joint(Y, l$lbar, l$l2bar, omega))
      dat$beta$beta_m <- beta_upd$m
      dat$beta$beta_s <- sqrt(diag(beta_upd$S))
      dat$beta$beta_var <- beta_upd$S
      dat$f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
    }
    if(dat$ll){
      myll <- c(myll, with(dat, calc_ell2(Y, l$lbar, l$l2bar, f$fbar, omega)))
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
  if(dat$ll){
    dat$myll <- myll
  }
  dat <- update_l_sequential(dat)
  dat$obj <- obj
  return(dat)
}
