ebmr_solve_future <- function(dat, max_iter, tol){

  fix_beta <- dat$fix_beta
  beta_joint <- dat$beta_joint
  #est_tau <- dat$est_tau

  check <- 1
  obj <- c()
  obj_old <- -Inf
  i <- 1

  #if(dat$est_tau){
  #  dat$omega_given <- dat$omega
  #}
  dat$obj_dec_warn <- FALSE
  while(i < max_iter & check > tol){
    #dat <- calc_lbar_future(dat)
    dat <- update_l_sequential_future(dat)
    #if(dat$ll){
    ll <- with(dat, calc_ell2(Y, l$lbar_o, l$l2bar_o, f$fbar, omega))
    obj <- c(obj, ll + dat$l$kl)
    #}
    # beta update
    if(!fix_beta & !beta_joint){
      dat <- update_beta_sequential_future(dat)
    }else if(!fix_beta & beta_joint){
      beta_upd <- with(dat,
                       update_beta_joint(Y, l$lbar, l$l2bar, omega))
      dat$beta$beta_m <- beta_upd$m
      dat$beta$beta_s <- sqrt(diag(beta_upd$S))
      dat$beta$beta_var <- beta_upd$S
      dat$f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
    }

    ll <- with(dat, calc_ell2(Y, l$lbar_o, l$l2bar_o, f$fbar, omega))

    obj <- c(obj, ll + dat$l$kl)

    obj_new <- obj[length(obj)]
    check <- obj_new - obj_old
    obj_old <- obj_new

    if(check < -1e-5){
      dat$obj_dec_warn <- TRUE
      #warning("Objective decreased, something is wrong.\n")
    }
    check <- abs(check)
    cat(i, ": ", obj_new, " ", dat$beta$beta_m, "\n")

    i <- i + 1
  }

  #dat <- update_l_sequential(dat)
  dat$obj <- obj

  return(dat)
}
