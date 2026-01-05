## Contains functions for fitting the ebnm model with point-normal prior with weights.


pn_precomp_weights <- function(x, s, w, par_init, fix_par) {
  fix_mu  <- fix_par[3]

  if (!fix_mu && any(s == 0)) {
    stop("The mode cannot be estimated if any SE is zero (the gradient does ",
         "not exist).")
  }

  if (any(s == 0)) {
    which_s0 <- which(s == 0)
    which_x_nz <- which(x[which_s0] != par_init$mu)
    i_x_nz <- which_s0[which_x_nz]
    #n0 <- length(which_s0) - length(which_x_nz)
    n0 <- sum(w[which_s0]) - sum(w[i_x_nz])
    #n1 <- length(which_x_nz)
    n1 <- sum(w[i_x_nz])
    sum1 <- sum(w[i_x_nz]*(x[i_x_nz] - par_init$mu)^2)
    x <- x[-which_s0]
    s <- s[-which_s0]
    w <- w[-which_s0]
  } else {
    n0 <- 0
    n1 <- 0
    sum1 <- 0
  }
  #n2 <- length(x)
  n2 <- sum(w)

  s2 <- s^2

  if (fix_mu) {
    z <- (x - par_init$mu)^2 / s2
    sum_z <- sum(w*z)
  } else {
    z <- NULL
    sum_z <- NULL
  }

  return(list(n0 = n0, n1 = n1, sum1 = sum1, n2 = n2, s2 = s2, z = z, sum_z = sum_z))
}

pn_nllik_weights <- function(par, x, s, w, par_init, fix_par,
                     n0, n1, sum1, n2, s2, z, sum_z,
                     calc_grad, calc_hess) {
  fix_pi0 <- fix_par[1]
  fix_s2  <- fix_par[2]
  fix_mu  <- fix_par[3]

  i <- 1
  if (fix_pi0) {
    alpha <- par_init$alpha
  } else {
    alpha <- par[i]
    i <- i + 1
  }
  if (fix_s2) {
    beta <- par_init$beta
  } else {
    beta <- par[i]
    i <- i + 1
  }
  if (fix_mu) {
    mu <- par_init$mu
  } else {
    mu <- par[i]
    z <- (x - mu)^2 / s2
    sum_z <- sum(w*z)
  }

  logist.alpha  <- 1 / (1 + exp(-alpha)) # scalar
  logist.nalpha <- 1 / (1 + exp(alpha))

  logist.beta   <- 1 / (1 + s2 * exp(-beta)) # scalar or vector
  logist.nbeta  <- 1 / (1 + exp(beta) / s2)

  y <- 0.5 * (z * logist.beta + log(logist.nbeta)) # vector

  # Negative log likelihood.
  C <- pmax(y, alpha)
  if (n0 == 0 || logist.alpha == 0) {
    nllik <- 0
  } else {
    nllik <- -n0 * log(logist.alpha)
  }
  nllik <- nllik - (n1 + n2) * (log(logist.nalpha))
  if (n1 > 0) {
    nllik <- nllik + 0.5 * n1 * beta
  }
  if (sum1 > 0) {
    nllik <- nllik + 0.5 * sum1 * exp(-beta)
  }
  nllik <- nllik + 0.5 * sum_z - sum(w*(log(exp(y - C) + exp(alpha - C)) + C))

  if (calc_grad || calc_hess) {
    dlogist.beta  <- logist.beta * logist.nbeta

    logist.y  <- 1 / (1 + exp(alpha - y)) # vector
    logist.ny <- 1 / (1 + exp(y - alpha))

    # Gradient.
    grad <- numeric(length(par))
    i <- 1
    if (!fix_pi0) {
      grad[i] <- -n0 * logist.nalpha + (n1 + n2) * logist.alpha - sum(w*logist.ny)
      i <- i + 1
    }
    if (!fix_s2) {
      dy.dbeta <- 0.5 * (z * dlogist.beta - logist.beta)
      grad[i] <- 0.5 * (n1 - sum1 * exp(-beta)) - sum(w*logist.y * dy.dbeta)
      i <- i + 1
    }
    if (!fix_mu) {
      dy.dmu <- (mu - x) * logist.beta / s2
      grad[i] <- sum(w*(mu - x) / s2) - sum(w*logist.y * dy.dmu)
    }
    attr(nllik, "gradient") <- grad
  }

  if (calc_hess) {
    dlogist.alpha <- logist.alpha * logist.nalpha
    dlogist.y <- logist.y * logist.ny

    # Hessian.
    hess <- matrix(nrow = length(par), ncol = length(par))
    i <- 1
    if (!fix_pi0) {
      tmp <- (n0 + n1 + n2) * dlogist.alpha
      hess[i, i] <- tmp - sum(w*dlogist.y)
      j <- i + 1
      if (!fix_s2) {
        hess[i, j] <- hess[j, i] <- sum(w*dlogist.y * dy.dbeta)
        j <- j + 1
      }
      if (!fix_mu) {
        hess[i, j] <- hess[j, i] <- sum(w*dlogist.y * dy.dmu)
      }
      i <- i + 1
    }
    if (!fix_s2) {
      d2y.dbeta2 <- 0.5 * ((z * (logist.nbeta - logist.beta) - 1) * dlogist.beta)
      tmp <- 0.5 * sum1 * exp(-beta) - sum(w*dlogist.y * dy.dbeta^2)
      hess[i, i] <- tmp - sum(w*logist.y * d2y.dbeta2)
      j <- i + 1
      if (!fix_mu) {
        d2y.dbetadmu <- (mu - x) * dlogist.beta / s2
        tmp <- -sum(w*dlogist.y * dy.dbeta * dy.dmu)
        hess[i, j] <- hess[j, i] <- tmp - sum(w*logist.y * d2y.dbetadmu)
      }
      i <- i + 1
    }
    if (!fix_mu) {
      tmp <- sum(1 / s2 - dlogist.y * dy.dmu^2)
      hess[i, i] <- tmp - sum(w*logist.y * logist.beta / s2)
    }
    attr(nllik, "hessian") <- hess
  }

  return(nllik)
}






mle_parametric_weights <- function(x,
                           s,
                           w,
                           par_init,
                           fix_par,
                           scalepar_fn,
                           precomp_fn,
                           nllik_fn,
                           postcomp_fn,
                           optmethod,
                           control,
                           use_grad,
                           use_hess) {
  #scale_factor <- 1 / median(s[s > 0])
  scale_factor <- 1
  x <- x * scale_factor
  s <- s * scale_factor

  par_init <- do.call(scalepar_fn, list(par = par_init,
                                        scale_factor = scale_factor))

  precomp <- do.call(precomp_fn, list(x = x,
                                      s = s,
                                      w = w,
                                      par_init = par_init,
                                      fix_par = fix_par))

  # Parameters that end up getting passed to all optimization functions:
  fn_params <- c(list(x = x, s = s, w = w, par_init = par_init, fix_par = fix_par),
                 precomp)

  p <- unlist(par_init)[!fix_par]

  # Don't initialize using infinite values.
  which.inf <- is.infinite(p)
  if (any(which.inf)) {
    p[which.inf] <- 0
  }

  if (all(fix_par)) {
    optpar <- par_init
    optval <- do.call(nllik_fn, c(list(par = NULL), fn_params,
                                  list(calc_grad = FALSE, calc_hess = FALSE)))
  } else if (optmethod == "nlm") {
    control <- modifyList(ebnm:::nlm_control_defaults(), control)

    optres <- do.call(nlm, c(list(f = nllik_fn, p = p),
                             fn_params,
                             list(calc_grad = use_grad, calc_hess = use_hess),
                             control))
    optpar <- optres$estimate
    optval <- optres$minimum
  } else if (optmethod == "trust") {
    control <- modifyList(trust_control_defaults(), control)

    # trust requires both a gradient and a Hessian.
    fn <- function(par, ...) {
      nllik <- do.call(nllik_fn,
                       list(par = par, calc_grad = TRUE, calc_hess = TRUE, ...))
      return(list(value = nllik,
                  gradient = attr(nllik, "gradient"),
                  hessian = attr(nllik, "hessian")))
    }
    optres <- do.call(trust::trust, c(list(objfun = fn, parinit = p),
                                      fn_params,
                                      control))
    optpar <- optres$argument
    optval <- optres$value
  } else if (optmethod == "lbfgsb") {
    control <- modifyList(lbfgsb_control_defaults(), control)

    # optim cannot accept a Hessian.
    fn <- function(par, ...) {
      return(do.call(nllik_fn,
                     list(par = par, calc_grad = FALSE, calc_hess = FALSE, ...)))
    }
    if (use_grad) {
      gr <- function(par, ...) {
        nllik <- do.call(nllik_fn,
                         list(par = par, calc_grad = TRUE, calc_hess = FALSE, ...))
        return(attr(nllik, "gradient"))
      }
    } else {
      gr <- NULL
    }

    optres <- do.call(optim, c(list(par = p, fn = fn, gr = gr),
                               fn_params,
                               list(control = control),
                               list(method = "L-BFGS-B")))
    optpar <- optres$par
    optval <- optres$value
  } else if (optmethod == "optimize") {
    control <- modifyList(optimize_control_defaults(), control)

    optres <- do.call(optimize, c(list(f = nllik_fn), fn_params,
                                  list(calc_grad = FALSE, calc_hess = FALSE),
                                  control))
    optpar <- optres$minimum
    optval <- optres$objective
  }

  # Combine the fixed and estimated parameters.
  retpar <- par_init
  retpar[!fix_par] <- optpar

  # Re-scale parameters and log likelihood.
  retpar <- do.call(scalepar_fn, list(par = retpar,
                                      scale_factor = 1 / scale_factor))
  optval <- optval - sum(w[is.finite(x)]) * log(scale_factor)
  #optval <- optval - sum(is.finite(x)) * log(scale_factor)

  retlist <- do.call(postcomp_fn, c(list(optpar = retpar,
                                         optval = optval,
                                         x = x,
                                         s = s,
                                         w = w,
                                         par_init = par_init,
                                         fix_par = fix_par,
                                         scale_factor = scale_factor),
                                    precomp))

  return(retlist)
}

# Postcomputations. A constant was subtracted from the log likelihood and needs
#   to be added back in. We also check boundary solutions here.
#
pn_postcomp_weights <- function(optpar, optval, x, s, w, par_init, fix_par, scale_factor,
                        n0, n1, sum1, n2, s2, z, sum_z) {
  llik <- pn_llik_from_optval_weights(optval, n1, n2, s2, w)
  retlist <- list(par = optpar, val = llik)

  # Check the solution pi0 = 1.
  fix_pi0 <- fix_par[1]
  fix_mu  <- fix_par[3]
  if (!fix_pi0 && fix_mu) {
    pi0_llik <- sum(w*-0.5 * log(2 * pi * s^2) - 0.5 * (w*(x - par_init$mu)^2 / s^2))
    pi0_llik <- pi0_llik + sum(w[is.finite(x)]) * log(scale_factor)
    if (pi0_llik > llik) {
      retlist$par$alpha <- Inf
      retlist$par$beta <- 0
      retlist$val <- pi0_llik
    }
  }

  return(retlist)
}


pn_llik_from_optval_weights <- function(optval, n1, n2, s2, w) {
  if (length(s2) == 1) {
    sum.log.s2 <- n2 * log(s2)
  } else {
    sum.log.s2 <- sum(w*log(s2))
  }

  return(-optval - 0.5 * ((n1 + n2) * log(2 * pi) + sum.log.s2))
}




parametric_workhorse_weights <- function(x,
                                 s,
                                 w,
                                 mode,
                                 scale,
                                 pointmass,
                                 g_init,
                                 fix_g,
                                 output,
                                 optmethod,
                                 control,
                                 checkg_fn,
                                 initpar_fn,
                                 scalepar_fn,
                                 precomp_fn,
                                 nllik_fn,
                                 postcomp_fn,
                                 summres_fn,
                                 partog_fn,
                                 postsamp_fn,
                                 call) {
  # I'm not sure why this is the case, but I run into infinite recursion issues
  #   if I don't extract the things I need from call here.
  call <- list(mode = call$mode, scale = call$scale)

  # Check that argument g_init is valid. All parametric families currently
  #   call into function check_g_init below.
  do.call(checkg_fn, list(g_init = g_init,
                          fix_g = fix_g,
                          mode = mode,
                          scale = scale,
                          pointmass = pointmass,
                          call = call))

  # Translate ebnm interface (mode/scale/pointmass/g_init/fix_g) into a
  #   generalized optimization interface (par_init/fix_par). fix_par, a vector
  #   of length 3, indicates whether 1. the weight of the spike component;
  #   2. the scale of the slab component; and 3. the location of the components
  #   is fixed.
  par_init <- do.call(initpar_fn, list(g_init = g_init,
                                       mode = mode,
                                       scale = scale,
                                       pointmass = pointmass,
                                       x = x,
                                       s = s))
  if (fix_g) {
    fix_par <- c(TRUE, TRUE, TRUE)
  } else {
    fix_par <- c(!pointmass,
                 !identical(scale, "estimate"),
                 !identical(mode, "estimate"))
  }

  optmethod <- ebnm:::handle_optmethod_parameter(optmethod, fix_par)

  # Don't use observations with infinite SEs when estimating g.
  x_optset <- x
  s_optset <- s
  w_optset <- w
  if (any(is.infinite(s))) {
    x_optset <- x[is.finite(s)]
    s_optset <- s[is.finite(s)]
    w_optset <- w[is.finite(s)]
  }

  # Estimate g. Function mle_parametric returns a list with fields par (which
  #   gives the estimated values of the parameters) and val (the optimal
  #   log likelihood attained).
  optres <- mle_parametric_weights(x = x_optset,
                           s = s_optset,
                           w = w_optset,
                           par_init = par_init,
                           fix_par = fix_par,
                           scalepar_fn = scalepar_fn,
                           precomp_fn = precomp_fn,
                           nllik_fn = nllik_fn,
                           postcomp_fn = postcomp_fn,
                           optmethod = optmethod$fn,
                           control = control,
                           use_grad = optmethod$use_grad,
                           use_hess = optmethod$use_hess)

  # Build return object.
  retlist <- list()

  if (ebnm:::data_in_output(output)) {
    retlist <- ebnm:::add_data_to_retlist(retlist, x, s)
  }

  if (ebnm:::posterior_in_output(output)) {
    posterior <- do.call(summres_fn, list(x = x,
                                          s = s,
                                          optpar = optres$par,
                                          output = output))
    retlist <- ebnm:::add_posterior_to_retlist(retlist, posterior, output, x)
  }

  if (ebnm:::g_in_output(output)) {
    fitted_g <- do.call(partog_fn, list(par = optres$par))
    retlist  <- ebnm:::add_g_to_retlist(retlist, fitted_g)
  }

  if (ebnm:::llik_in_output(output)) {
    loglik  <- optres$val
    retlist <- ebnm:::add_llik_to_retlist(retlist, loglik, x, df = sum(!fix_par))
  }

  if (ebnm:::sampler_in_output(output)) {
    sampler <- function(nsamp) {
      samp <- postsamp_fn(x, s, optres$par, nsamp)
      colnames(samp) <- names(x)
      return(samp)
    }
    retlist <- ebnm:::add_sampler_to_retlist(retlist, sampler)
  }

  return(retlist)
}

ebnm_pn_weights <- function(x, w, s = 1, mode = 0,
                            scale = "estimate",
                            g_init = NULL,
                            fix_g = FALSE,
                            output = ebnm_output_default(),
                            optmethod = NULL,
                            control = NULL,
                            ...){
  call <- match.call()
  if(is.null(control)){
    control <- list()
  }
  retlist <- parametric_workhorse_weights(x = x,
                                  s = s,
                                  w = w,
                                  mode = mode,
                                  scale = scale,
                                  pointmass = TRUE,
                                  g_init = g_init,
                                  fix_g = fix_g,
                                  output = output,
                                  optmethod = optmethod,
                                  control = control,
                                  checkg_fn = ebnm:::pn_checkg,
                                  initpar_fn = ebnm:::pn_initpar,
                                  scalepar_fn = ebnm:::pn_scalepar,
                                  precomp_fn = pn_precomp_weights,
                                  nllik_fn = pn_nllik_weights,
                                  postcomp_fn = pn_postcomp_weights,
                                  summres_fn = ebnm:::pn_summres,
                                  partog_fn = ebnm:::pn_partog,
                                  postsamp_fn = ebnm:::pn_postsamp,
                                  call = call,
                                  ...)
  return(ebnm:::as_ebnm(retlist, call))
}


############ Test

# pn_nll_slow <- function(par, x, s, w){
#
#   alpha <- par[1]
#   beta <- par[2]
#   mu <- par[3]
#
#   logist.alpha  <- 1 / (1 + exp(-alpha)) # scalar
#   logist.nalpha <- 1 / (1 + exp(alpha))
#
#   sigma2 <- exp(beta)
#
#   log_lik <- logist.alpha*dnorm(x = x, mean = 0, sd = s) +
#              logist.nalpha*dnorm(x = x, mean = mu, sd = sqrt(s^2 + sigma2))
#   log_lik <- log(log_lik)
#   log_lik <- sum(w*log_lik)
#   return(log_lik)
# }
#
# pn_nll_ebnm <- function(par, x, s, w, fix_par = c(FALSE, FALSE, TRUE), use_scale = TRUE){
#
#   if(use_scale){
#     scale_factor <- 1 / median(s[s > 0])
#   }else{
#     scale_factor <- 1
#   }
#   x <- x * scale_factor
#   s <- s * scale_factor
#
#   par_init <- list(alpha = par[1],
#                    beta = par[2],
#                    mu = par[3])
#
#   par_init <- do.call(ebnm:::pn_scalepar, list(par = par_init,
#                                                scale_factor = scale_factor))
#
#   par <- unlist(par_init[!fix_par])
#   precomp <- pn_precomp_weights(x, s, w, par_init, fix_par= fix_par)
#
#   fn_params <- c(list(par = par, x = x, s = s, w = w,
#                       par_init = par_init, fix_par = fix_par, calc_grad = TRUE, calc_hess = TRUE),
#                  precomp)
#
#   ll <- do.call(pn_nllik_weights, fn_params)
#
#   optval <- ll[1]
#   retpar <- par_init
#   retpar[!fix_par] <- par
#
#
#   # Re-scale parameters and log likelihood.
#   retpar <- do.call(ebnm:::pn_scalepar, list(par = retpar,
#                                              scale_factor = 1 / scale_factor))
#
#   optval <- optval - sum(w[is.finite(x)]) * log(scale_factor)
#
#   retlist <- do.call(pn_postcomp_weights, c(list(optpar = retpar,
#                                                  optval = optval,
#                                                  x = x,
#                                                  s = s,
#                                                  w = w,
#                                                  par_init = par_init,
#                                                  fix_par = fix_par,
#                                                  scale_factor = scale_factor),
#                                                  precomp))
#
#   return(retlist)
# }
#
#
#
#
# set.seed(2026)
# n <- 1e4
# z <- rbinom(n = n, size = 1, prob = 0.3)
# theta <- rep(0, n)
# theta[z == 1] <- rnorm(n = sum(z == 1), mean = 0, sd = 5)
#
# s <- runif(n = n, min = 0.5, max = 5)
#
# x <- rnorm(n = n, mean = theta, sd = s)
# w <- runif(n = n)
#
# a <- log(0.2/0.8)
# pn_nll_slow(par = c(a, log(8), 0), x = x, s =s, w = w)
#
# pn_nll_ebnm(par = c(a, log(8), 0), x = x, s =s, w = w, use_scale= FALSE)
# pn_nll_ebnm(par = c(a, log(8), 0), x = x, s =s, w = w, use_scale= TRUE)
