
#'@title Empirical Shrinkage Multivariable MR
#'@param beta_hat_Y Vector of SNP-outcome associations (length p)
#'@param se_Y Standard errors of beta_hat_Y
#'@param beta_hat_X Matrix of SNP-exposure associations (p by K)
#'@param se_X matrix of standard errors of beta_hat_X
#'@param G G matrix. If NULL, G will be estimated using the method given in g_type.
#'@param R Optional correlation matrix for overlapping samples.
#'@param ebnm_fn Options prior distribution family. Defaults to point-normal.
#'@param max_iter Maximum number of iterations
#'@param sigma_beta Optional prior variance for causal parameters
#'@param tol Convergence tolerance
#'@param pval_thresh p-value threshold (suggest 1)
#'@param beta_joint Use joint updates for beta (suggest TRUE)
#'@param g_type Method to estimate G. Suggest "gfa"
#'@param augment_G Augment estimated G
#'@export
esmr <- function(beta_hat_Y, se_Y,
                 beta_hat_X, se_X,
                 G = NULL,
                 R = NULL,
                 ebnm_fn = flashier:::as.ebnm.fn(prior_family = "point_normal", optmethod = "nlm"),
                 max_iter = 100,
                 sigma_beta = Inf,
                 tol = default_precision(c(ncol(beta_hat_X)+1, nrow(beta_hat_X))),
                 pval_thresh =1,
                 nrand = 0,
                 #post_prob_thresh = 0,
                 lfsr_thresh = 0.01,
                 beta_m_init = NULL,
                 which_beta = NULL,
                 fix_beta = FALSE,
                 beta_joint = TRUE,
                 g_type = c("gfa", "svd"),
                 svd_zthresh = 0,
                 augment_G = FALSE){


  #if(length(fix_beta) > 1 & beta_joint) stop("if beta_joint = TRUE, fix_beta should have length 1.\n")
  g_type <- match.arg(g_type, choices = c("gfa", "svd"))
  dat <- set_data(beta_hat_Y, se_Y, beta_hat_X, se_X, R)


  if(is.null(G)){
    if(dat$p == 2){
      G <- diag(dat$p)
    }else if(!missing(which_beta)){
      warning("Cannot estimate G for network problem yet.\n")
      G <- diag(dat$p)
    }else{
      G <- estimate_G(beta_hat_X <- dat$Y[,-1],
                      se_X = dat$S[,-1],
                      R = R[-1, -1, drop = FALSE],
                      type = g_type,
                      svd_zthresh = svd_zthresh,
                      augment = augment_G)
    }
  }
  dat$G <- check_matrix(G, "G", n = dat$p)
  dat$k <- ncol(dat$G)

  dat$beta <- init_beta(dat$p, which_beta, beta_m_init, fix_beta)
  dat$f <- make_f(dat)


  dat$pval_thresh <- pval_thresh
  #dat$post_prob_thresh <- post_prob_thresh

  if(pval_thresh < 1){
    pval <- with(dat, 2*pnorm(-abs(Y/S)))
    pval_min <- apply(pval[,-1,drop = FALSE], 1, min)
    ix <- which(pval_min < pval_thresh)
    if(nrand > 0){
      nrand <- min(nrand, dat$n - length(ix))
      nullix <- sample(setdiff(1:dat$n, ix), size = nrand, replace = FALSE)
      ix <- sort(c(ix, nullix))
    }
    dat <- subset_data(dat, ix)
  }

  dat$l <- init_l(dat$n, dat$k, ncol(dat$G))


  dat$beta_joint <- beta_joint
  dat$ebnm_fn <- ebnm_fn
  dat$sigma_beta <- sigma_beta
  dat$lfsr_thresh <- lfsr_thresh

  #dat$est_tau <- est_tau

  dat <- esmr_solve(dat, max_iter, tol )

  return(dat)
}




