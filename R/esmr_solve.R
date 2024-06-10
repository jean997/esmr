esmr_solve <- function(dat, max_iter, tol){

  check <- 1
  obj <- c()
  obj_old <- -Inf
  i <- 1

  dat$obj_dec_warn <- FALSE
  nb <- length(dat$beta$beta_j)
  ixlist <- list()
  while(i < max_iter & check > tol){
    # l update
    dat <- update_l_sequential(dat, seq(dat$p), dat$g_init, dat$fix_g)
    dat <- update_l_sequential(dat, seq(dat$p), dat$g_init, dat$fix_g)

    ll <- with(dat, calc_ell2(Y, l$abar, l$a2bar, f$fgbar, omega))
    obj <- c(obj, ll + dat$l$kl)

    # beta update
    if(!dat$beta_joint){
      dat <- update_beta_sequential(dat)
    }else{
      dat$beta$V <- matrix(0, nrow = nb, ncol = nb)
      jj <- unique(dat$beta$beta_j)
      for(j in jj){
        ii <- which(dat$beta$beta_j == j & !dat$beta$fix_beta)
        if(length(ii) == 0) next
        ix <- dat$beta$beta_k[ii]
        beta_upd <- update_beta_joint(dat, j = j, ix = ix)
        #beta_upd <- update_beta_joint(du, j = j, ix = ix) ## temporary

        dat$beta$beta_m[ii] <- beta_upd$m
        dat$beta$beta_s[ii] <- sqrt(diag(beta_upd$S))
        dat$beta$V[ii,ii] <- beta_upd$S
        dat$f <- make_f(dat)
      }
    }
    ll <- with(dat, calc_ell2(Y, l$abar, l$a2bar, f$fgbar, omega))
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
  dat$ixlist <- ixlist
  return(dat)
}
