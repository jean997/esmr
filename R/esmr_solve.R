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
    dat <- update_l_sequential(dat)
    dat <- update_l_sequential(dat)

    ll <- with(dat, calc_ell2(Y, l$abar, l$a2bar, f$fgbar, omega))
    obj <- c(obj, ll + dat$l$kl)

    # beta update
    if(!dat$beta_joint){
      dat <- update_beta_sequential(dat)
      dat$beta$V <- diag(dat$beta$beta_s^2)
    }else{
      dat$beta$V <- matrix(0, nrow = nb, ncol = nb)
      if(dat$R_is_id | length(unique(dat$beta$beta_j)) == 1){
        # if all omega are diagonal or only estimating one row, update F by rows
        jj <- unique(dat$beta$beta_j)
        for(j in jj){
          ii <- which(dat$beta$beta_j == j & !dat$beta$fix_beta)
          if(length(ii) == 0) next
          ix <- dat$beta$beta_k[ii]
          beta_upd <- update_beta_joint(dat, j = j, ix = ix)

          dat$beta$beta_m[ii] <- beta_upd$m
          dat$beta$beta_s[ii] <- sqrt(diag(beta_upd$S))
          dat$beta$V[ii,ii] <- beta_upd$S
        }
      }else{
        e_ix <- which(!dat$beta$fix_beta)
        ub <- update_beta_full_joint(dat, prior_cov = NULL)
        dat$beta$beta_m[e_ix] <- ub$m
        dat$beta$V[e_ix,e_ix] <- ub$S
        dat$beta$beta_s[e_ix] <- sqrt(diag(ub$S))
      }
      dat$f <- make_f(dat)
    }
    ## new step, update total effects based on constraints
    if(!is.null(dat$which_tot_c)){
      f <- t(complete_T(t(dat$f$fbar), dat$which_tot_c)$total_effects)
      ## temporary fix for negative variance
      f2 <- pmax(t(complete_T(t(dat$f$f2bar), dat$which_tot_c)$total_effects), f^2) ## assumes independent estimates
      vf <- f2 - (f^2)
      ix <- cbind(dat$beta$beta_j, dat$beta$beta_k)
      dat$beta$beta_m <- f[ix]
      dat$beta$beta_s <- sqrt((f2[ix]) - (f[ix])^2)
      diag(dat$beta$V) <- dat$beta$beta_s^2
      dat$f <- make_f(dat)
    }

    ###
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
