
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
                 ebnm_fn = flashier::flash_ebnm(prior_family = "point_normal", optmethod = "nlm"),
                 max_iter = 100,
                 sigma_beta = Inf,
                 tol = default_precision(c(ncol(beta_hat_X)+1, nrow(beta_hat_X))),
                 #####
                 direct_effect_template = NULL,
                 direct_effect_init = NULL,
                 # add ability to fix some effects later
                 #fix_beta = FALSE,
                 beta_joint = TRUE,
                 g_type = c("gfa", "svd"),
                 svd_zthresh = 0,
                 augment_G = TRUE,
                 ix1 = NULL){


  #if(length(fix_beta) > 1 & beta_joint) stop("if beta_joint = TRUE, fix_beta should have length 1.\n")
  g_type <- match.arg(g_type, choices = c("gfa", "svd"))
  stopifnot(beta_joint %in% c(TRUE, FALSE))

  dat <- set_data(beta_hat_Y, se_Y, beta_hat_X, se_X, R)
  dat$beta_joint <- beta_joint
  dat$ebnm_fn <- ebnm_fn
  dat$sigma_beta <- sigma_beta



  if(is.null(G)){
    if(dat$p == 2){
      G <- diag(dat$p)
    }else if(!missing(which_beta)){
      warning("Cannot estimate G for network problem yet.\n")
      G <- diag(dat$p)
    }else{
      G <- estimate_G(beta_hat_X <- dat$Y[,-1,drop =F],
                      se_X = dat$S[,-1, drop = F],
                      R = R[-1, -1, drop = FALSE],
                      type = g_type,
                      svd_zthresh = svd_zthresh,
                      augment = augment_G)
    }
  }
  dat$G <- check_matrix(G, "G", n = dat$p)

  dat$k <- ncol(dat$G)


  if(!is.null(direct_effect_template)){
    B <- check_B_template(direct_effect_template)
    which_beta <- rbind(B$which_tot_u, B$which_tot_c)[,c(2,1)] ## transpose
    colnames(which_beta) <- c("row", "col")
    if(!is.null(direct_effect_init)){ ## to do: add check for valid initialization
      beta_m_init <- direct_effect_init(which_beta)
    }else{
      beta_m_init <- NULL
    }
    fix_beta <- c(rep(FALSE, nrow(B$which_tot_u)), rep(TRUE, nrow(B$which_tot_c)))
    if(nrow(B$which_tot_c) > 0){
      dat$which_tot_c <- B$which_tot_c # no transpose[,c(2, 1)]
      #colnames(dat$which_tot_c) <- c("row", "col")
    }
  }else{
    which_beta <- beta_m_init <- fix_beta <- NULL
  }

  dat$beta <- init_beta(dat$p, which_beta, beta_m_init, fix_beta)
  dat$f <- make_f(dat)
  #dat$f0 <- make_f(dat)

  dat$l <- init_l(dat$n, dat$p, dat$k)

  if(is.null(ix1)){
    dat <- esmr_solve(dat, max_iter, tol )
  }else{
    dat <- get_ix1_ix0(dat, ix1)
    dat <- subset_data(dat, dat$ix1)
    dat <- esmr_solve(dat, max_iter, tol)
  }

  return(dat)
}




