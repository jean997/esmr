esmr_resolve <- function(
    mod, direct_effect_template = NULL, max_iter = 100,  tol = "default",
    restrict_dag = TRUE) {

  if(tol == "default"){
    tol <- default_precision(c(ncol(dat$Y), nrow(dat$Y)))
  }

  mod$B_template <- direct_effect_template
#  browser()
  B <- check_B_template(mod$B_template, mod$p, restrict_dag = restrict_dag)
  #
  which_beta <- rbind(B$which_tot_u, B$which_tot_c)[,c(2,1), drop=FALSE] ## transpose
  which_beta <- cbind(which_beta, c(rep(FALSE, nrow(B$which_tot_u)), rep(TRUE, nrow(B$which_tot_c))))
  colnames(which_beta) <- c("row", "col", "fixed")

  # Only keep betas that are in the new template
  original_idx <- data.frame(row = mod$beta$beta_j, col = mod$beta$beta_k, idx = seq_along(mod$beta$beta_j))
  keep_betas <- merge(original_idx, which_beta, by = c("col", "row"))
  keep_betas <- arrange(keep_betas, fixed)

#  keep_betas <- mod$beta$beta_j == which_beta[,1] & mod$beta$beta_k == which_beta[,2]
  mod$beta$beta_j <- mod$beta$beta_j[keep_betas$idx]
  mod$beta$beta_k <- mod$beta$beta_k[keep_betas$idx]
  mod$beta$beta_m <- mod$beta$beta_m[keep_betas$idx]
  mod$beta$beta_s <- mod$beta$beta_s[keep_betas$idx]
  mod$beta$V <- mod$beta$V[keep_betas$idx, keep_betas$idx]
  mod$beta$fix_beta <- keep_betas$fixed

  mod <- esmr_solve(mod, max_iter, tol)

  ### Pasted from esmr
  # TODO: Could refactor into separate function
  ## post-processing
  #o <- match(1:mod$p, mod$traits)
  #mod <- reorder_mod(mod, o)

  if (!is.null(direct_effect_template) && restrict_dag) {
    mod$direct_effects <- total_to_direct(t(mod$f$fbar) - diag(mod$p))
    delt_pvals <- delta_method_pvals(mod)
    mod$pvals_dm <- delt_pvals$pmat
    mod$se_dm <- delt_pvals$semat
  }

  # Reformat beta_hat and beta_se to matrix format
  beta_hat <- beta_se <- matrix(0, nrow = mod$p, ncol = mod$p)
  fix_beta <- matrix(FALSE, nrow = mod$p, ncol = mod$p)
  # Lower triangular format
  beta_ind <- cbind(mod$beta$beta_k, mod$beta$beta_j)
  beta_hat[beta_ind] <- mod$beta$beta_m
  # TODO: Need to fix this since beta_se is not the same as beta_m
  beta_se[beta_ind] <- mod$beta$beta_s
  fix_beta[beta_ind] <- mod$beta$fix_beta

  mod$beta_mat <- list(
    beta_hat = beta_hat,
    beta_se = beta_se
  )

  mod$elbo <- tail(mod$obj, n = 1)

  return(mod)
}
