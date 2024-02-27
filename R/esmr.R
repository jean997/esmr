
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
esmr <- function(beta_hat_X, se_X,
                 beta_hat_Y=NULL, se_Y = NULL,
                 G = NULL,
                 R = NULL,
                 pval_thresh = 5e-8,
                 variant_ix = NULL,
                 ###
                 ebnm_fn = flashier::flash_ebnm(prior_family = "point_normal", optmethod = "nlm"),
                 max_iter = 100,
                 sigma_beta = Inf,
                 tol = "default",
                 #####
                 direct_effect_template = NULL,
                 direct_effect_init = NULL,
                 # add ability to fix some effects later
                 # direct_effect_fix = NULL,
                 #fix_beta = FALSE,
                 beta_joint = TRUE,
                 augment_G = TRUE){


  #if(length(fix_beta) > 1 & beta_joint) stop("if beta_joint = TRUE, fix_beta should have length 1.\n")
  #g_type <- match.arg(g_type, choices = c("gfa", "svd"))
  g_type <- "gfa"
  dat <- set_data(beta_hat_Y, se_Y, beta_hat_X, se_X, R)
  stopifnot(beta_joint %in% c(TRUE, FALSE))

  dat <- set_data(beta_hat_Y, se_Y, beta_hat_X, se_X, R)
  if(is.null(G)){
    if(dat$p == 2){
      G <- diag(dat$p)
    }else if(!missing(direct_effect_template) & !is.null(direct_effect_template)){
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
  dat$G <- check_matrix(G, n = dat$p)
  dat <- order_upper_tri(dat, direct_effect_template, direct_effect_init)
  dat <- init_beta(dat)
  dat$beta_joint <- beta_joint
  dat$ebnm_fn <- ebnm_fn
  dat$sigma_beta <- sigma_beta
  dat$R_is_id <- is.null(R) | all(R == diag(dat$p))

  dat$k <- ncol(dat$G)

  dat$f <- make_f(dat)

  dat$l <- init_l(dat$n, dat$p, dat$k)



  dat$beta_joint <- beta_joint
  dat$ebnm_fn <- ebnm_fn
  dat$sigma_beta <- sigma_beta
  #dat$lfsr_thresh <- lfsr_thresh

  if(!is.null(pval_thresh)){
    dat <- get_ix1_ix0(dat, paste0("pval-", pval_thresh))
    dat <- subset_data(dat, dat$ix1)
  }else if(!is.null(variant_ix)){
    dat <- subset_data(dat, variant_ix)
  }
  if(tol == "default"){
    tol <- default_precision(c(ncol(dat$Y), nrow(dat$Y)))
  }
  dat <- esmr_solve(dat, max_iter, tol)


  o <- match(1:dat$p, dat$traits)
  dat <- reorder_data(dat, o)

  dat$direct_effects <- total_to_direct(t(dat$f$fbar) - diag(dat$p))
  delt_pvals <- delta_method_pvals(dat)
  dat$pvals_dm <- delt_pvals$pmat
  dat$se_dm <- delt_pvals$semat
  return(dat)
}




