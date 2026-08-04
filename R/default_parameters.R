mvmr_default_params <- function(){
  list(
    tau_init = NULL,
    fix_tau = FALSE,
    ebnm_fn = flashier::flash_ebnm(prior_family = "point_normal", optmethod = "nlm"),
    g_init = NULL,
    fix_g = FALSE,
    tol = "default",
    beta_prior_cov = NULL,
    beta_joint = TRUE,
    augment_G = TRUE,
    cond_num = 1e10,
    strict_mode = TRUE,
    keep_ebnm_res = FALSE
  )
}


nesmr_default_params <- function(){
  list(
    direct_effect_init = NULL,
    restrict_dag = TRUE,
    tau_init = NULL,
    fix_tau = FALSE,
    ebnm_fn = flashier::flash_ebnm(prior_family = "point_normal", optmethod = "nlm"),
    g_init = NULL,
    fix_g = FALSE,
    tol = "default",
    beta_prior_cov = NULL,
    beta_joint = TRUE,
    cond_num = 1e10,
    strict_mode = TRUE,
    keep_ebnm_res = FALSE
  )
}

factor_default_params <- function(){
  list(
    tau_init = NULL,
    fix_tau = FALSE,
    ebnm_fn = flashier::flash_ebnm(prior_family = "point_normal", optmethod = "nlm"),
    g_init = NULL,
    fix_g = FALSE,
    tol = "default",
    beta_prior_cov = NULL,
    beta_joint = TRUE,
    cond_num = 1e10,
    keep_ebnm_res = FALSE
  )
}

