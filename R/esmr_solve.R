esmr_solve <- function(dat, max_iter, tol){

  check <- 1
  obj <-  c()
  obj_old <- -Inf
  i <- 1

  dat$obj_dec_warn <- FALSE
  cond_num <- dat$cond_num
  if(is.null(cond_num)) cond_num <- 1e10

  nb <- length(dat$beta$beta_j)
  low_info_flag <- FALSE
  while(i < max_iter && check > tol && !low_info_flag){
    low_info_flag <- FALSE

    # l update
    dat <- update_l_sequential(dat, seq(dat$k), dat$g_init, dat$fix_g)
    #dat <- update_l_sequential(dat, seq(dat$p), dat$g_init, dat$fix_g)

    ll <- with(dat, calc_ell2(Y, l$abar, l$a2bar, f$fgbar, omega, omega_logdet, s_equal))
    obj <- c(obj, ll + dat$l$kl + dat$beta$kl)

    # beta update
    remove_suggest <- NULL
    if(!dat$beta_joint){
      dat <- update_beta_sequential(dat)
      dat$beta$V <- diag(dat$beta$beta_s^2)
    }else{
      dat$beta$V <- matrix(0, nrow = nb, ncol = nb)
      if((dat$R_is_id | length(unique(dat$beta$beta_j)) == 1) & !dat$is_factors){
        # if all omega are diagonal or only estimating one row, update F by rows
        jj <- unique(dat$beta$beta_j)
        for(j in jj){
          ii <- which(dat$beta$beta_j == j & !dat$beta$fix_beta)
          if(length(ii) == 0) next
          ix <- dat$beta$beta_k[ii]
          beta_upd <- update_beta_joint(dat, j = j, ix = ix, ii = ii, cond_num = cond_num)

          dat$beta$beta_m[ii] <- beta_upd$m
          dat$beta$beta_s[ii] <- sqrt(diag(beta_upd$S))
          dat$beta$V[ii,ii] <- beta_upd$S
          if(dat$is_factors){
            dat$f <- make_f_factors(dat)
          }else{
            dat$f <- make_f(dat)
          }
          if(!is.null(beta_upd$remove_suggest)) {
            low_info_flag <- TRUE
            remove_suggest <- beta_upd$remove_suggest
          }
        }
      }else{
        e_ix <- which(!dat$beta$fix_beta)
        ub <- update_beta_full_joint(dat, cond_num = cond_num)
        dat$beta$beta_m[e_ix] <- ub$m
        dat$beta$V[e_ix,e_ix] <- ub$S
        dat$beta$beta_s[e_ix] <- sqrt(diag(ub$S))
        if(dat$is_factors){
          dat$f <- make_f_factors(dat)
        }else{
          dat$f <- make_f(dat)
        }
        if(!is.null(ub$remove_suggest)) {
          low_info_flag <- TRUE
          remove_suggest <- ub$remove_suggest
        }
      }
    }

    # Update KL divergence for beta if we have a prior
    if(!is.null(dat$beta$prior_cov) && length(dat$beta$prior_cov) > 0){
      kl_ix <- !dat$beta$fix_beta
      prior_cov_mat <- dat$beta$prior_cov * diag(sum(kl_ix))
      # Note: Can pass prior_precision instead to avoid solving a bunch of times
      dat$beta$kl <- - kl_mvn(
        dat$beta$beta_m[kl_ix], dat$beta$V[kl_ix, kl_ix,drop=F], 0, prior_cov_mat)
    }else{
      dat$beta$kl <- 0
    }

    ## update total effects based on constraints
    if(any(dat$beta$fix_beta)){
      which_const <- cbind(dat$beta$beta_k, dat$beta$beta_j)[dat$beta$fix_beta,,drop = FALSE]
      colnames(which_const) <- c("row", "col")
      f <- t(complete_T(t(dat$f$fbar), which_const)$total_effects)
      ix <- cbind(dat$beta$beta_j, dat$beta$beta_k)
      dat$beta$beta_m <- f[ix]

      if(dat$is_factors){
        dat$f <- make_f_factors(dat)
      }else{
        dat$f <- make_f(dat)
      }
    }

    ## tau update
    if(!is.null(dat$tau) & !dat$fix_tau){
      min_tau <- dat$tau/10
      max_tau <- dat$tau*10
      if(dat$tau == 0){
        max_tau <- 10*median(dat$S^2)
      }
      dat <- update_tau(dat,tau_min = min_tau, tau_max = max_tau)
      #ll <- with(dat, calc_ell2(Y, l$abar, l$a2bar, f$fgbar, omega, omega_logdet, s_equal))
      obj <- c(obj, ll + dat$l$kl + dat$beta$kl)
    }

    ###
    ll <- with(dat, calc_ell2(Y, l$abar, l$a2bar, f$fgbar, omega, omega_logdet, s_equal))
    #cat("ll: ", ll, "l$kl: ", dat$l$kl, "beta$kl: ", dat$beta$kl, "\n")
    obj <- c(obj, ll + dat$l$kl + dat$beta$kl)

    obj_new <- obj[length(obj)]
    check <- obj_new - obj_old
    #check <- max(abs(dat$beta$beta_m - beta_old))
    obj_old <- obj_new
    #beta_old <- dat$beta$beta_m
    #cat(check, "\n")
    if(check < -1e-12){
      dat$obj_dec_warn <- TRUE
      warning("Objective decreased, something may be wrong.\n")
    }
    cat(i, ": ", obj_new, " ", dat$beta$beta_m, " ", dat$tau, "\n")
    #cat(i, ": ", check, " ", dat$beta$beta_m, "\n")

    i <- i + 1

    if(dat$is_factors && low_info_flag){
      dat <- remove_worst_factor(dat)
    } else if (low_info_flag) {
      dat$remove_suggest <- remove_suggest
    }
  }

  dat$obj <- obj

  return(dat)
}
