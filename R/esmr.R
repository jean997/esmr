
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
                 pval_thresh = 5e-8,
                 variant_ix = NULL,
                 ebnm_fn = flashier::flash_ebnm(prior_family = "point_normal", optmethod = "nlm"),
                 g_init = NULL,
                 fix_g = FALSE,
                 max_iter = 100,
                 sigma_beta = Inf,
                 tol = "default",
                 #####
                 beta_m_init = NULL,
                 which_beta = NULL,
                 fix_beta = FALSE,
                 beta_joint = TRUE,
                 augment_G = TRUE){


  #if(length(fix_beta) > 1 & beta_joint) stop("if beta_joint = TRUE, fix_beta should have length 1.\n")
  #g_type <- match.arg(g_type, choices = c("gfa", "svd"))
  g_type <- "gfa"
  dat <- set_data(beta_hat_Y, se_Y, beta_hat_X, se_X, R)
  stopifnot(beta_joint %in% c(TRUE, FALSE))

  if(is.null(G)){
    if(dat$p == 2){
      G <- diag(dat$p)
    }else if(!missing(which_beta)){
      warning("Cannot estimate G for network problem yet.\n")
      G <- diag(dat$p)
    }else{
      G <- estimate_G(beta_hat_X = dat$Y[,-1,drop =F],
                      se_X = dat$S[,-1, drop = F],
                      R = R[-1, -1, drop = FALSE],
                      type = g_type,
                      augment = augment_G)
    }
  }
  dat$G <- check_matrix(G, "G", n = dat$p)

  dat$k <- ncol(dat$G)

  dat$beta <- init_beta(dat$p, which_beta, beta_m_init, fix_beta)
  dat$f <- make_f(dat)
  #dat$f0 <- make_f(dat)

  dat$l <- init_l(dat$n, dat$p, dat$k)


  dat$beta_joint <- beta_joint
  dat$ebnm_fn <- ebnm_fn
  dat$g_init <- g_init
  dat$fix_g <- fix_g

  dat$sigma_beta <- sigma_beta
  #dat$lfsr_thresh <- lfsr_thresh

  if(!is.null(variant_ix)){
    dat <- subset_data(dat, variant_ix)
  }else if(!is.null(pval_thresh)){
    dat <- get_ix1_ix0(dat, paste0("pval-", pval_thresh))
    dat <- subset_data(dat, dat$ix1)
  }
  if(tol == "default"){
    tol <- default_precision(c(ncol(dat$Y), nrow(dat$Y)))
  }
  dat <- esmr_solve(dat, max_iter, tol)

  return(dat)
}




