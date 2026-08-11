
#'@title ESMR Workhorse Called from wrapper functions esmr, nesmr, and esmr_factors
#'@export
esmr_workhorse <- function(beta_hat_X, se_X,
                 beta_hat_Y=NULL, se_Y = NULL,
                 beta_hat_Z=NULL, se_Z = NULL,
                 factors_matrix = NULL,
                 factors_residual_sd = NULL,
                 G = NULL,
                 R = NULL,
                 pval_thresh = NULL,
                 variant_ix = NULL,
                 ld_scores = NULL,
                 RE = NULL,
                 tau_init = NULL,
                 fix_tau = FALSE,
                 ###
                 ebnm_fn = flashier::flash_ebnm(prior_family = "point_normal", optmethod = "trust"),
                 g_init = NULL,
                 fix_g = FALSE,
                 max_iter = 100,
                 tol = "default",
                 restrict_dag = TRUE,
                 #####
                 direct_effect_template = NULL,
                 direct_effect_init = NULL,
                 # add ability to fix some effects later
                 # direct_effect_fix = NULL,
                 #fix_beta = FALSE,
                 strict_mode = FALSE,
                 beta_prior_cov = NULL,
                 beta_joint = TRUE,
                 total_to_direct_bootstrap_samples = 10000,
                 keep_ebnm_res = FALSE,
                 augment_G = TRUE,
                 cond_num = 1e10){


  #if(length(fix_beta) > 1 & beta_joint) stop("if beta_joint = TRUE, fix_beta should have length 1.\n")
  #g_type <- match.arg(g_type, choices = c("gfa", "svd"))
  g_type <- "gfa"
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

  # note set_data returns object of class c("esmr", "list")
  if(!is.null(beta_hat_Z)){
    dat <- set_data_factors(beta_hat_Y, se_Y, beta_hat_X, se_X,
                            beta_hat_Z, se_Z, factors_matrix, factors_residual_sd,
                            R, ld_scores, RE, tau_init)
    dat$is_factors <- TRUE
  }else{
    dat <- set_data(beta_hat_Y, se_Y, beta_hat_X, se_X, R, ld_scores, RE, tau_init)
    dat$is_factors <- FALSE
  }

  dat$is_nesmr <- !is.null(direct_effect_template)
  if(dat$is_nesmr){
    class(dat) <- c("nesmr", class(dat))
  }

  dat$R_is_id <- (is.null(R) || all(R == diag(dat$p))) & is.null(RE)
  dat$cond_num <- cond_num
  dat$beta_joint <- beta_joint
  dat$ebnm_fn <- ebnm_fn
  dat$g_init <- g_init
  dat$fix_g <- fix_g
  dat$fix_tau <- fix_tau
  dat$restrict_dag <- restrict_dag
  dat$strict_mode <- strict_mode

  if (dat$is_nesmr) {
    class(dat) <- c(c("nesmr"), class(dat))
  }

  if(is.null(G) & !dat$is_nesmr & !dat$is_factors){
    if(dat$p == 2){
      dat$G <- diag(dat$p)
    } else{
      dat$G <- estimate_G(beta_hat_X = dat$Y[,-1,drop =F],
                      se_X = dat$S[,-1, drop = F],
                      R = R[-1, -1, drop = FALSE],
                      type = g_type,
                      augment = augment_G)
    }
  }else if(dat$is_nesmr){
    dat$G <- diag(1, dat$p)
  }else if(dat$is_factors){
    dat$G <- diag(1, dat$k)
  }else{
    dat$G <- check_matrix(G, n = dat$p)
  }
  dat$k <- ncol(dat$G)

  if(dat$is_nesmr){
    dat <- order_upper_tri(dat, direct_effect_template, direct_effect_init,
                         restrict_dag = restrict_dag)
  }

  dat <- init_beta(dat)
  if (!is.null(beta_prior_cov)) {
    nb <- sum(! dat$beta$fix_beta)
    dat$beta$prior_cov <- check_beta_prior_cov(beta_prior_cov, nb)
    if (length(dat$beta$prior_cov) > 0) {
        dat$beta$prior_precision <- solve(dat$beta$prior_cov)
    }

  }

  if(dat$is_factors){
    dat$f <- make_f_factors(dat)
    dat$l <- init_l(dat$n, dat$k, dat$k)
  }else{
    dat$f <- make_f(dat)
    dat$l <- init_l(dat$n, dat$p, dat$k)
  }

  # subset variants
  if(!is.null(variant_ix)){
    dat <- subset_data(dat, variant_ix)
  }else if(!is.null(pval_thresh)){
    dat <- get_ix1_ix0(
      dat,
      paste0("pval-", pval_thresh),
      remove_empty_B_cols = dat$is_nesmr)

    dat <- subset_data(dat, dat$ix1)
  }
  if(tol == "default"){
    tol <- default_precision(c(ncol(dat$Y), nrow(dat$Y)))
  }

  ## solve esmr problem
  dat <- esmr_solve(dat, max_iter, tol, keep_ebnm_res = keep_ebnm_res)

  if (strict_mode && !is.null(dat$remove_suggest)) {
    stop(sprintf("Strict mode is on and a low information trait was suggested for removal. Consider removing trait %s and re-running ESMR.", dat$remove_suggest))
  }

  ## post-processing
  if(dat$is_nesmr){
    o <- match(1:dat$p, dat$traits)
    dat <- reorder_data(dat, o)
  }

  if (dat$is_nesmr && is_dag(dat$B_template) && !all(direct_effect_template == 0)) {
    # Multiply by direct effect template to ensure rounding is not an issue
    dat$direct_effects <- total_to_direct(t(dat$f$fbar) - diag(dat$p)) * dat$B_template
    delt_pvals <- delta_method_pvals(dat)
    dat$direct_effects_se <- delt_pvals$semat * dat$B_template
    dat$direct_effects_log_pval <- delt_pvals$pmat * dat$B_template
    dat$standard_error_method <- "delta_method"
  } else if (dat$is_nesmr && !is_dag(dat$B_template)) {
    dat$total_effects <- t(dat$f$fbar)
    dat$direct_effects <- total_to_direct(dat$total_effects - diag(dat$p), restrict_dag = FALSE) * dat$B_template

    bootstrap_ix <- cbind(
      dat$beta$beta_k,
      dat$beta$beta_j
    )
    pbstrap <- total_to_direct_parameteric_bootstrap(
      x = dat$beta$beta_m,
      beta_cov = dat$beta$V,
      beta_ix = bootstrap_ix,
      d = ncol(dat$total_effects),
      bootstrap_samples = 1000
    )

    dat$direct_effect_se <- pbstrap$direct_effects_se
    dat$direct_effects_bootstrap <- pbstrap$direct_effects_bootstrap
    dat$standard_error_method <- "parametric_bootstrap"
  }

  if (dat$is_factors) {
    fbar <- dat$f$fbar
    nfactor <- dat$k - 2
    total_effect_matrix <- cbind(t(fbar[1:2, , drop = FALSE]),
                                 rbind(matrix(0, nrow = 2, ncol = nfactor), diag(nfactor)))
    dat$direct_effects <- t(total_to_direct(total_effect_matrix - diag(dat$k)))
  }
  dat <- format_betas(dat)
  dat$elbo <- tail(dat$obj, n = 1)

  return(dat)
}




