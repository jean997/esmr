ebmr_solve_future <- function(dat, max_iter, tol){

  check <- 1
  obj <- c()
  obj_old <- -Inf
  i <- 1

  dat$obj_dec_warn <- FALSE
  nb <- length(dat$beta$beta_j)
  while(i < max_iter & check > tol){
    # l update
    dat <- update_l_sequential_future(dat)

    ll <- with(dat, calc_ell2(Y, l$lbar_o, l$l2bar_o, f$fbar, omega))
    obj <- c(obj, ll + dat$l$kl)

    # beta update
    if(!dat$beta_joint){
      dat <- update_beta_sequential_future(dat)
    }else{
      V <- matrix(0, nrow = nb, ncol = nb)
      jj <- unique(dat$beta$beta_j)
      for(j in jj){
        ii <- which(dat$beta$beta_j == j & !dat$beta$fix_beta)
        if(length(ii) == 0) next
        ix <- dat$beta$beta_k[ii]
        beta_upd <- update_beta_joint_future(dat, j = j, ix = ix)
        dat$beta$beta_m[ii] <- beta_upd$m
        dat$beta$beta_s[ii] <- sqrt(diag(beta_upd$S))
        V[ii,ii] <- beta_upd$S

        dat$f <- make_f(dat)
      }
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

  dat$obj <- obj

  return(dat)
}
