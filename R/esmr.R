
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
esmr <- function(beta_hat_X, se_X,
                 beta_hat_Y=NULL, se_Y = NULL,
                 G = NULL,
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
                 restrict_dag = TRUE,
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
  stopifnot(beta_joint %in% c(TRUE, FALSE))
  if(! is.null(pval_thresh) && ! is.null(variant_ix)){
    stop("Please specify only one of pval_thresh or variant_ix.")
  }

  is_nesmr <- !is.null(direct_effect_template)
  if(!is.null(ld_scores) | !is.null(RE)){
    if(is.null(ld_scores) | is.null(RE)){
      stop("Please specify both ld_scores and RE to include correction for GWAS confounding.")
    }
    if(is.null(tau_init)){
      tau_init <- 1e-4
    }
  }

  dat <- set_data(beta_hat_Y, se_Y, beta_hat_X, se_X, R, ld_scores, RE, tau_init)

  class(dat) <- c(c("esmr"), class(dat))
  dat$is_nesmr <- ! is.null(direct_effect_template)
  if (dat$is_nesmr) {
    class(dat) <- c(c("nesmr"), class(dat))
  }

  if(is.null(G)){
    if(dat$p == 2){
      G <- diag(dat$p)
    } else if (dat$is_nesmr) {
      warning("Cannot estimate G for network problem yet.\n")
      G <- diag(dat$p)
    } else{
      G <- estimate_G(beta_hat_X = dat$Y[,-1,drop =F],
                      se_X = dat$S[,-1, drop = F],
                      R = R[-1, -1, drop = FALSE],
                      type = g_type,
                      augment = augment_G)
    }
  }
  dat$G <- check_matrix(G, n = dat$p)

  dat <- order_upper_tri(dat, direct_effect_template, direct_effect_init,
                         restrict_dag = restrict_dag)

  dat <- init_beta(dat, restrict_dag = restrict_dag)
  dat$beta_joint <- beta_joint
  dat$ebnm_fn <- ebnm_fn
  dat$sigma_beta <- sigma_beta
  dat$R_is_id <- (is.null(R) || all(R == diag(dat$p))) & is.null(RE)

  dat$k <- ncol(dat$G)

  dat$f <- make_f(dat)

  dat$l <- init_l(dat$n, dat$p, dat$k)

  dat$g_init <- g_init
  dat$fix_g <- fix_g

  dat$fix_tau <- fix_tau
  dat$restrict_dag <- restrict_dag

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

  ## post-processing
  o <- match(1:dat$p, dat$traits)
  dat <- reorder_data(dat, o)

  if (!is.null(direct_effect_template) && restrict_dag) {
    dat$direct_effects <- total_to_direct(t(dat$f$fbar) - diag(dat$p))
    delt_pvals <- delta_method_pvals(dat)
    dat$pvals_dm <- delt_pvals$pmat
    dat$se_dm <- delt_pvals$semat
  }

  # Reformat beta_hat and beta_se to matrix format
  beta_hat <- beta_se <- matrix(0, nrow = dat$p, ncol = dat$p)
  fix_beta <- matrix(FALSE, nrow = dat$p, ncol = dat$p)
  # Lower triangular format
  beta_ind <- cbind(dat$beta$beta_k, dat$beta$beta_j)
  beta_hat[beta_ind] <- dat$beta$beta_m
  beta_se[beta_ind] <- dat$beta$beta_s
  fix_beta[beta_ind] <- dat$beta$fix_beta

  dat$beta_mat <- list(
    beta_hat = beta_hat,
    beta_se = beta_se
  )

  dat$elbo <- tail(dat$obj, n = 1)

  return(dat)
}




