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

      # if (dat$logdet_penalty) {
      #   B_tot <- matrix(0, nrow = dat$p, ncol = dat$p)
      #   ix <- cbind(dat$beta$beta_j, dat$beta$beta_k)
      #   e_ix <- which(!dat$beta$fix_beta)
      #   B_tot[ix] <- dat$beta$beta_m
      #
      #   result <- optim(
      #     par = as.vector(B_tot),
      #     fn = h_det_vector,
      #     gr = h_det_grad_vector,
      #     s = 1,
      #     d = dat$p,
      #     method = "BFGS",
      #     control = list(maxit = 100),
      #     hessian = TRUE
      #   )
      #
      #   if (result$convergence != 0) {
      #     new_f_mat <- matrix(result$par, nrow = dat$p, ncol = dat$p)
      #     # Set beta to the new values
      #     dat$beta$beta_m <- new_f_mat[ix[e_ix, ]]
      #   } else {
      #     browser()
      #     warning("Optimization did not converge.")
      #   }

        # TODO: Not sure what to do about the standard errors...

      dat$f <- make_f(dat)
    }

    # TODO: How do we update the standard deviations...?
    # I think we can do this here...

    ## new step, update total effects based on constraints
    if(any(dat$beta$fix_beta)){
      which_const <- cbind(dat$beta$beta_k, dat$beta$beta_j)[dat$beta$fix_beta,,drop = FALSE]
      colnames(which_const) <- c("row", "col")
      f <- t(complete_T(t(dat$f$fbar), which_const)$total_effects)
      ## assumes independent estimates
      f2_1 <- t(complete_T(t(dat$f$f2bar), which_const)$total_effects)
      f2_2 <- f^2
      f2_3 <- t(complete_T(t(dat$f$fbar)^2, which_const)$total_effects)
      f2 <- f2_1 + f2_2 - f2_3 ## E[f^2] = g(E[f^2*]) + g(f^2*) - g(E[f*]^2)
      ###
      vf <- f2 - (f^2)
      ix <- cbind(dat$beta$beta_j, dat$beta$beta_k)
      dat$beta$beta_m <- f[ix]
      dat$beta$beta_s <- sqrt((f2[ix]) - (f[ix])^2)
      diag(dat$beta$V) <- dat$beta$beta_s^2
      dat$f <- make_f(dat)
    }

    # ## Minimize the objective function with h_det penalty
    # if (dat$logdet_penalty) {
    #   # B_tot <- matrix(0, nrow = dat$p, ncol = dat$p)
    #   ix <- cbind(dat$beta$beta_j, dat$beta$beta_k)
    #   #e_ix <- which(!dat$beta$fix_beta)
    #   # B_tot[ix] <- dat$beta$beta_m
    #
    #   #         make_f
    #
    #
    #   # TODO: Write function that is - ell2 with betas used for fgbar
    #   # Keep everything else the same..
    #   # optim_func <- function(x) {
    #   #   calc_ell2(Y, l$abar, dat$l$a2bar, f$fgbar, omega)
    #   # }
    #
    #   ell_func <- function(fg) {
    #     X <- matrix(fg, ncol = dat$p)
    #     - calc_ell2(dat$Y, dat$l$abar, dat$l$a2bar, X, dat$omega)
    #   }
    #
    #   result <- optim(
    #     par = as.vector(dat$f$fgbar - diag(dat$p)),
    #     fn = h_det_ell,
    #     gr = h_det_ell_grad,
    #     s = 1,
    #     d = dat$p,
    #     ell_func = ell_func,
    #     method = "BFGS",
    #     control = list(maxit = 100),
    #     #hessian = TRUE
    #   )
    #
    #   if (result$convergence != 0) {
    #     new_f_mat <- matrix(result$par, nrow = dat$p, ncol = dat$p)
    #     # Set beta to the new values
    #     dat$beta$beta_m <- new_f_mat[ix]
    #   } else {
    #     browser()
    #     warning("Optimization did not converge.")
    #   }
    # }

    dat$f <- make_f(dat)

    ###
    # Could we use "calc_ell2" as the objective function to maximize?
    ll <- with(dat, calc_ell2(Y, l$abar, l$a2bar, f$fgbar, omega))
    obj <- c(obj, ll + dat$l$kl)

    obj_new <- obj[length(obj)]
    check <- obj_new - obj_old
    obj_old <- obj_new

    if(check < -1e-5){
      dat$obj_dec_warn <- TRUE
      warning("Objective decreased, something is wrong.\n")
      next
    }
    check <- abs(check)
    cat(i, ": ", obj_new, " ", dat$beta$beta_m, "\n")

    i <- i + 1
  }

  dat$obj <- obj
  dat$ixlist <- ixlist
  return(dat)
}
