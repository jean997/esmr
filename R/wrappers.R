#'@title Empirical Shrinkage Multivariable MR
#'@param beta_hat_Y Vector of SNP-outcome associations (length n)
#'@param se_Y Standard errors of beta_hat_Y
#'@param beta_hat_X Matrix of SNP-exposure associations (n by p-1).
#'@param se_X matrix of standard errors of beta_hat_X (n by p-1).
#'@param G G matrix. If NULL, G will be estimated.
#'@param R Nuisance correlation matrix if there is sample overlap (optional).
#'@param pval_thresh p-value threshold for variant selection (optional).
#'@param variant_ix Indices of variants to include (optional).
#'@param ld_scores List of ld_scores. If missing, no additional variance term will be added.
#'@param RE RE
#'@param max_iter Maximum number of iterations.
#'@param params Additional parameters. See mvmr_default_params() for more information. Generally these can be left at their defaults.
#'@export
esmr <- function(beta_hat_X, se_X,
                 beta_hat_Y=NULL, se_Y = NULL,
                 G = NULL,
                 R = NULL,
                 pval_thresh = NULL,
                 variant_ix = NULL,
                 ld_scores = NULL,
                 RE = NULL,
                 max_iter = 100,
                 params = mvmr_default_params()){

  args <- as.list(environment())
  args["params"] <- NULL
  default_params <- mvmr_default_params()

  for(n in names(default_params)){
    if(is.null(params[[n]])) params[[n]] <- default_params[[n]]
  }
  for(n in names(params)){
    if(! n %in% names(default_params)) stop("Unknown parameter ", n, " provided.")
  }
  args <- c(args, params)
  call_esmr_workhorse(args)
}

#'@title Network Empirical Shrinkage MR
#'@param beta_hat_X Matrix of SNP-trait associations (n by p)
#'@param se_X Matrix of standard errors of beta_hat_X
#'@param direct_effect_template p by p matrix with all elements equal to 0 or 1 specifying the DAG structure.
#'@param R Nuisance correlation matrix if there is sample overlap (optional).
#'@param pval_thresh p-value threshold for variant selection (optional).
#'@param variant_ix Indices of variants to include (optional).
#'@param ld_scores List of ld_scores. If missing, no additional variance term will be added.
#'@param RE RE
#'@param max_iter Maximum number of iterations.
#'@param params Additional parameters. See nesmr_default_params() for more information. Generally these can be left at their defaults.
#'@export
nesmr <- function(beta_hat_X, se_X,
                  direct_effect_template,
                  R = NULL,
                  pval_thresh = NULL,
                  variant_ix = NULL,
                  ld_scores = NULL,
                  RE = NULL,
                  max_iter = 100,
                  params = nesmr_default_params()){
  args <- as.list(environment())
  args["params"] <- NULL
  default_params <- nesmr_default_params()

  for(n in names(default_params)){
    if(is.null(params[[n]])) params[[n]] <- default_params[[n]]
  }
  for(n in names(params)){
    if(! n %in% names(default_params)) stop("Unknown parameter ", n, " provided.")
  }
  args <- c(args, params)
  call_esmr_workhorse(args)
}

#'@title Empirical Shrinkage Multivariable MR with Factors
#'@param beta_hat_Y Vector of SNP-outcome associations (length n)
#'@param se_Y Standard errors of beta_hat_Y
#'@param beta_hat_X Vector of SNP-exposure associations (length n)
#'@param se_X Standard errors of beta_hat_X
#'@param beta_hat_Z Matrix of SNP-confounder associations (n by p-2).
#'@param se_Z matrix of standard errors of beta_hat_Z (n by p-2).
#'@param factors_matrix Factors matrix (p-2 by k).
#'@param factors_residual_sd Factors residual sd
#'@param R Nuisance correlation matrix if there is sample overlap (optional).
#'@param pval_thresh p-value threshold for variant selection (optional).
#'@param variant_ix Indices of variants to include (optional).
#'@param ld_scores List of ld_scores. If missing, no additional variance term will be added.
#'@param RE RE
#'@param max_iter Maximum number of iterations.
#'@param params Additional parameters. See factor_default_params() for more information. Generally these can be left at their defaults.
#'@export
esmr_factors <- function(beta_hat_X, se_X,
                         beta_hat_Y, se_Y,
                         beta_hat_Z, se_Z,
                         factors_matrix,
                         factors_residual_sd = NULL,
                         R = NULL,
                         pval_thresh = NULL,
                         variant_ix = NULL,
                         ld_scores = NULL,
                         RE = NULL,
                         max_iter = 100,
                         params = factor_default_params()){
  args <- as.list(environment())
  args["params"] <- NULL
  default_params <- factor_default_params()

  for(n in names(default_params)){
    if(is.null(params[[n]])) params[[n]] <- default_params[[n]]
  }
  for(n in names(params)){
    if(! n %in% names(default_params)) stop("Unknown parameter ", n, " provided.")
  }
  args <- c(args, params)
  call_esmr_workhorse(args)
}

call_esmr_workhorse <- function(args) {
  rlang::try_fetch(
    rlang::inject(esmr_workhorse(!!!args)),
    error = function(cnd) {
      rlang::abort(
        message = "",#conditionMessage(cnd),
        # This points the error specifically to the workhorse call
        call = rlang::call2("esmr_workhorse", !!!args),
        # Chaining allows rlang::last_trace() to see the relationship
        parent = cnd,
        # This hides the 'call_esmr_workhorse' internal plumbing
        arg = "args"
      )
    }
  )
}