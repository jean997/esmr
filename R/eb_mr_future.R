
#'@title Empirical Shrinkage Multivariable MR
#'@param beta_hat_Y Vector of SNP-outcome associations (length p)
#'@param se_Y Standard errors of beta_hat_Y
#'@param beta_hat_X Matrix of SNP-exposure associations (p by K)
#'@param se_X matrix of standard errors of beta_hat_X
#'@param R Optional correlation matrix for overlapping samples.
#'@param ebnm_fn Options prior distribution family. Defaults to point-normal.
#'@param pval_thresh p-value threshold for beta update
#'@export
eb_mr_future <- function(beta_hat_Y, se_Y,
                         beta_hat_X, se_X,
                         G = NULL,
                  R = NULL,
                  ebnm_fn = flashier:::as.ebnm.fn(prior_family = "point_normal", optmethod = "nlm"),
                  max_iter = 100,
                  sigma_beta = Inf,
                  tol = default_precision(c(ncol(beta_hat_X)+1, nrow(beta_hat_X))),
                  pval_thresh =1,
                  post_prob_thresh = 0,
                  beta_m_init = NULL,
                  which_beta = NULL,
                  fix_beta = FALSE,
                  beta_joint = TRUE,
                  g_type = c("gfa", "svd"),
                  svd_zthresh = 0){
                  #est_tau = FALSE,
                  #ll = FALSE){

  g_type <- match.arg(g_type, g_type)
  dat <- set_data(beta_hat_Y, se_Y, beta_hat_X, se_X, R)

  dat$pval_thresh <- pval_thresh
  dat$post_prob_thresh <- post_prob_thresh

  if(pval_thresh < 1){
    pval <- with(dat, 2*pnorm(-abs(Y/S)))
    pval_min <- apply(pval[,-1,drop = FALSE], 1, min)
    ix <- which(pval_min < pval_thresh)
    dat <- subset_data(dat, ix)
  }else if(post_prob_thresh > 0){
    wpost <- get_wpost(dat$Y, dat$S, 2:p, prior_family = "point_normal")
    wpost_max <- apply(wpost, 1, max)
    ix <- which(wpost_max >= post_prob_thresh)
    dat <- subset_data(dat, ix)
  }

  n <- dat$n
  p <- dat$p
  if(is.null(G)){
    if(p == 2){
      G <- diag(p)
    }else if(!missing(which_beta)){
      warning("Cannot estimate G for network problem yet.\n")
      G <- diag(p)
    }else{
      G <- estimate_G(beta_hat_X <- dat$Y[,-1],
                      se_X = dat$S[,-1],
                      R = R,
                      type = g_type,
                      svd_zthresh = svd_zthresh)
    }
  }
  dat$G <- check_matrix(G, "G", n = p)
  k <- ncol(dat$G)

  dat$beta <- init_beta(p, which_beta, beta_joint, beta_m_init)
  dat$f_fun <- make_f_fun_future(p, dat$beta$beta_j, dat$beta$beta_k, G)

  dat$l <- init_l_future(n, k)
  dat$f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
  dat$fix_beta <- fix_beta
  dat$beta_joint <- beta_joint
  dat$ebnm_fn <- ebnm_fn
  dat$sigma_beta <- sigma_beta

  #dat$est_tau <- est_tau
  #dat$ll <- ll

  dat <- ebmr_solve_future(dat, max_iter, tol )

  return(dat)
}




