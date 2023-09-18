esmr_solve_2part <- function(dat, max_iter, tol){

  check <- 1
  obj <- c()
  obj_old <- -Inf
  i <- 1

  dat$obj_dec_warn <- FALSE
  nb <- length(dat$beta$beta_j)
  ixlist <- list()
  while(i < max_iter & check > tol){
    # l update
    #dat <- update_l_sequential(dat)
    dat <- update_l_sequential_2parts(dat)


    d1 <- subset_data(dat, dat$ix1)
    d0 <- subset_data(dat, dat$ix0)

    ll1 <- with(d1, calc_ell2(Y, l$abar, l$a2bar, f$fgbar, omega))
    ll0 <- with(d0, calc_ell2(Y, l$abar, l$a2bar, f0$fgbar, omega))
    obj <- c(obj, ll1 + ll0 + dat$l$kl)

    ### Temporary Code
    # minlfsr <- apply(dat$l$lfsr, 1, min)
    # ixk <- which(minlfsr < dat$lfsr_thresh)
    # du <- subset_data(dat, ixk)
    # ixlist[[i]] <- ixk
    ###

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
        #beta_upd <- update_beta_joint(dat, j = j, ix = ix)
        #beta_upd <- update_beta_joint(du, j = j, ix = ix) ## temporary
        beta_upd1 <- update_beta_joint(d1, j = j, ix = ix)
        beta_upd0 <- update_beta_joint(d0, j = j, ix = ix)

        dat$beta$beta_m[ii] <- beta_upd1$m
        dat$beta$beta_s[ii] <- sqrt(diag(beta_upd1$S))
        dat$beta$V[ii,ii] <- beta_upd1$S

        d0$beta$beta_m[ii] <- beta_upd0$m
        d0$beta$beta_s[ii] <- sqrt(diag(beta_upd0$S))

        dat$f0 <- make_f(d0)
        dat$f <- make_f(dat)
      }
    }
    d1 <- subset_data(dat, dat$ix1)
    d0 <- subset_data(dat, dat$ix0)
    ll1 <- with(d1, calc_ell2(Y, l$abar, l$a2bar, f$fgbar, omega))
    ll0 <- with(d0, calc_ell2(Y, l$abar, l$a2bar, f0$fgbar, omega))
    obj <- c(obj, ll1 + ll0 + dat$l$kl)
    #ll <- with(dat, calc_ell2(Y, l$abar, l$a2bar, f$fgbar, omega))

    #obj <- c(obj, ll + dat$l$kl)
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
  #dat$ixlist <- ixlist
  return(dat)
}
