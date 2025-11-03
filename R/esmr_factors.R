
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
#'@param pval_thresh p-value threshold for estimation
#'@param variant_ix Instead of using pval_thresh, directly specify the indices of variants used for estimation.
#'@param beta_joint Use joint updates for beta (suggest TRUE)
#'@param g_type Method to estimate G. Suggest "gfa"
#'@param augment_G Augment estimated G
#'@export
esmr_factors <- function(beta_hat_X, se_X,
                 beta_hat_Y=NULL, se_Y = NULL,
                 beta_hat_Z=NULL, se_Z = NULL,
                 factors_matrix = NULL,
                 factors_residual_sd = NULL,
                 R = NULL,
                 pval_thresh = NULL,
                 variant_ix = NULL,
                 ld_scores = NULL,
                 RE = NULL,
                 tau_init = NULL,
                 fix_tau = FALSE,
                 ###
                 ebnm_fn = flashier::flash_ebnm(prior_family = "point_normal", optmethod = "nlm"),
                 g_init = NULL,
                 fix_g = FALSE,
                 max_iter = 100,
                 sigma_beta = Inf,
                 tol = "default",
                 cond_num = 1e10,
                 ####
                 beta_joint = TRUE){



  stopifnot(beta_joint %in% c(TRUE, FALSE))
  if(! is.null(pval_thresh) && ! is.null(variant_ix)){
    stop("Please specify only one of pval_thresh or variant_ix.")
  }

  if(!is.null(ld_scores) | !is.null(RE)){
    if(is.null(ld_scores) | is.null(RE)){
      stop("Please specify both ld_scores and RE to include correction for GWAS confounding.")
    }
    if(is.null(tau_init)){
      tau_init <- 1e-4
    }
  }

  dat <- set_data_factors(beta_hat_Y, se_Y, beta_hat_X, se_X,
                  beta_hat_Z, se_Z, factors_matrix, factos_residual_sd,
                  R, ld_scores, RE, tau_init)

  dat$G <- diag(dat$k )
  dat$cond_num <- cond_num
  class(dat) <- c(c("esmr"), class(dat))
  dat$is_nesmr <- FALSE
  dat$is_factors <- TRUE

  dat <- init_beta_factors(dat)
  dat$beta_joint <- beta_joint
  dat$ebnm_fn <- ebnm_fn
  dat$sigma_beta <- sigma_beta
  dat$R_is_id <- (is.null(R) || all(R == diag(dat$p))) & is.null(RE)

  dat$f <- make_f_factors(dat)

  dat$l <- init_l(dat$n, dat$k , dat$k)

  dat$g_init <- g_init
  dat$fix_g <- fix_g

  dat$fix_tau <- fix_tau

  # subset variants
  if(!is.null(variant_ix)){
    dat <- subset_data(dat, variant_ix)
  }else if(!is.null(pval_thresh)){
    dat <- get_ix1_ix0(
      dat,
      paste0("pval-", pval_thresh),
      remove_empty_B_cols = is_nesmr)

    dat <- subset_data(dat, dat$ix1)
  }
  if(tol == "default"){
    tol <- default_precision(c(ncol(dat$Y), nrow(dat$Y)))
  }

  # Pre-compute log(det(omega))
  dat$omega_logdet <- get_omega_logdet(dat$omega, dat$s_equal, n = dat$n)

  ## solve esmr problem
  dat <- esmr_solve(dat, max_iter, tol)

  dat <- format_betas(dat)


  dat$elbo <- tail(dat$obj, n = 1)

  return(dat)
}




