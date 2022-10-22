
#'@title Empirical Shrinkage Multivariable MR
#'@param beta_hat_Y Vector of SNP-outcome associations (length p)
#'@param se_Y Standard errors of beta_hat_Y
#'@param beta_hat_X Matrix of SNP-exposure associations (p by K)
#'@param se_X matrix of standard errors of beta_hat_X
#'@param R Optional correlation matrix for overlapping samples.
#'@param ebnm_fn Options prior distribution family. Defaults to point-normal.
#'@param pval_thresh p-value threshold for beta update
#'@export
eb_mr <- function(beta_hat_Y, se_Y, beta_hat_X, se_X,
                  R = NULL,
                  ebnm_fn = flashier:::as.ebnm.fn(prior_family = "point_normal", optmethod = "nlm"),
                  max_iter = 100,
                  seed = 123, sigma_beta = Inf,
                  tol = default_precision(c(ncol(beta_hat_X)+1, nrow(beta_hat_X))),
                  pval_thresh =1, lfsr_thresh = 1,
                  beta_m_init = NULL, which_beta = NULL,
                  fix_beta = FALSE,
                  beta_joint = TRUE,
                  est_tau = FALSE,
                  ll = FALSE){

  set.seed(seed)

  dat <- set_data(beta_hat_Y, se_Y, beta_hat_X, se_X, R)
  n <- dat$n
  p <- dat$p

  dat$beta <- init_beta(p, which_beta, beta_joint, beta_m_init)
  dat$f_fun <- make_f_fun(p, dat$beta$beta_j, dat$beta$beta_k)

  dat$l <- init_l(n, p)
  dat$f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
  dat$fix_beta <- fix_beta
  dat$beta_joint <- beta_joint
  dat$ebnm_fn <- ebnm_fn
  dat$sigma_beta <- sigma_beta

  dat$pval <- with(dat, 2*pnorm(-abs(Y/S)))
  dat$pval_min <- apply(dat$pval[,-1,drop = FALSE], 1, min)
  dat$ix <- which(dat$pval_min < pval_thresh)

  dat$lfsr_thresh <- lfsr_thresh
  dat$pval_thresh <- pval_thresh
  dat$est_tau <- est_tau
  dat$ll <- ll

  dat <- ebmr_solve(dat, max_iter, tol )

  return(dat)
}


calc_ell <- function(Y, lbar, l2bar, fbar, f2bar, omega){
  n <- nrow(Y)
  p <- ncol(Y)
  stopifnot(nrow(lbar) == n & ncol(lbar) == p)
  stopifnot(nrow(l2bar) == n & ncol(l2bar) == p)
  stopifnot(nrow(fbar) == p & nrow(f2bar) == p & ncol(fbar) == p & ncol(f2bar) == p)
  #stopifnot(length(beta_m) == p-1 & length(beta_s) == p-1)


  if("matrix" %in% class(omega)){
    s_equal <- TRUE
  }else{
    stopifnot(class(omega) == "list")
    stopifnot(length(omega) == n)
    s_equal <- FALSE
  }

  ybar <- lbar %*% t(fbar)

  varlbar <- l2bar - (lbar^2)
  varfbar <- f2bar - (fbar^2)
  if(s_equal){
    part_a <- map_dbl(seq(nrow(ybar)), function(i){
      crossprod(t(Y[i,,drop =F]), omega) %>% tcrossprod(ybar[i, , drop = F])
    }) %>% sum()
    part_b <- map_dbl(seq(nrow(ybar)), function(i){
      crossprod(t(ybar[i,,drop =F]), omega) %>% tcrossprod(ybar[i, , drop = F])
    }) %>% sum()
    x1 <- sum(t( t(tcrossprod(lbar^2 , varfbar))*diag(omega)))
    x2 <- sum(t( t(tcrossprod(varlbar , fbar^2))*diag(omega)))
    x3 <- sum(t( t(tcrossprod(varlbar , varfbar))*diag(omega)))
    part_b <- part_b + x1 + x2 + x3
    ell <- part_a - 0.5*part_b
  }else{
    part_a <- map_dbl(seq(nrow(ybar)), function(i){
      crossprod(t(Y[i,,drop =F]), omega[[i]]) %>% tcrossprod(ybar[i, , drop = F])
    }) %>% sum()
    part_b <- map_dbl(seq(nrow(ybar)), function(i){
      crossprod(t(ybar[i,,drop =F]), omega[[i]]) %>% tcrossprod(ybar[i, , drop = F])
    }) %>% sum()
    x1 <- map_dbl(seq(nrow(ybar)), function(i){
      sum(t(t( tcrossprod((lbar^2)[i,,drop =F], varfbar))*diag(omega[[i]])) )
    }) %>% sum()
    x2 <- map_dbl(seq(nrow(ybar)), function(i){
      sum(t(t( tcrossprod(varlbar[i,,drop =F], fbar^2))*diag(omega[[i]])) )
    }) %>% sum()
    x3 <- map_dbl(seq(nrow(ybar)), function(i){
     sum(t(t( tcrossprod(varlbar[i,,drop =F], varfbar))*diag(omega[[i]])) )
    }) %>% sum()
    part_b <- part_b + x1 + x2 + x3
    ell <- part_a - 0.5*part_b
  }
  return(ell)
}

