default_precision <- function(dims){
  sqrt(.Machine$double.eps)*prod(dims)
}

make_f <- function(dat){
  if(length(dat$beta$beta_j) == 0 | is.null(dat$beta$beta_j)){
    return(list(fbar = diag(dat$p), #f2bar = diag(dat$p),
                fgbar = dat$G)) #, fg2bar = dat$G^2))
  }
  fbar <-  diag(dat$p)
  nb <- length(dat$beta$beta_j)
  ix <- cbind(dat$beta$beta_j, dat$beta$beta_k)

  fbar[ix] <- dat$beta$beta_m
  fgbar <- fbar %*% dat$G

  # V <- matrix(0, nrow = dat$p, ncol = dat$p)
  # V[ix] <- dat$beta$beta_s^2
  # f2bar <- (fbar^2) + V
  # fg2bar <- (fgbar^2) + (V %*% (dat$G^2))
  return(list(fgbar = fgbar, #fg2bar = fg2bar,
              fbar = fbar))#,f2bar = f2bar))
}

make_f_factors <- function(dat){

  fbar <- matrix(0, nrow = dat$p, ncol = dat$k )
  fbar[1,1] <- fbar[2,2] <- 1
  fbar[3:dat$p, 3:(dat$k )] <- dat$factors_matrix

  nb <- length(dat$beta$beta_j)
  ix <- cbind(dat$beta$beta_j, dat$beta$beta_k)

  fbar[ix] <- dat$beta$beta_m
  fgbar <- fbar

  return(list(fgbar = fgbar, #fg2bar = fg2bar,
              fbar = fbar))#,f2bar = f2bar))
}

format_betas <- function(dat){
  # Reformat beta_hat and beta_se to matrix format
  beta_hat <- beta_se <- matrix(0, nrow = dat$p, ncol = dat$k)
  fix_beta <- matrix(FALSE, nrow = dat$p, ncol = dat$k)
  # Lower triangular format
  beta_ind <- cbind(dat$beta$beta_j, dat$beta$beta_k)
  beta_hat[beta_ind] <- dat$beta$beta_m
  beta_se[beta_ind] <- dat$beta$beta_s
  fix_beta[beta_ind] <- dat$beta$fix_beta

  dat$beta$beta_hat <- beta_hat
  dat$beta$beta_se <- beta_se
  dat$beta$beta_fixed <- fix_beta
  return(dat)
}




get_omega <- function(R, S, s_equal, any_missing){

  p <- ncol(S)
  R_is_id <- is.null(R) | all(R == diag(p))

  if(s_equal && R_is_id){
    s <- S[1,]
    omega <- diag(1/s^2)
  }else if(s_equal){
    s <- S[1,]
    omega <- solve_diag_psd_diag(R, s)
  }else if(R_is_id && !any_missing){
    omega <- apply(S, 1, function(s){
      diag(1/s^2, nrow = p)
    }, simplify = FALSE)
  }else if(!any_missing){
    omega <- apply(S, 1, function(s){
      solve_diag_psd_diag(R, s)
    }, simplify = FALSE)
  }else if(R_is_id){
    S[is.na(S)] <- Inf
    omega <- apply(S, 1, function(s){
      diag(1/s^2, nrow = p)
    }, simplify = FALSE)
  }else{
    pat <- data.frame(is.na(S))
    pat$pat <- apply(pat, 1, function(x){paste0(x, collapse = "-")})
    pat_sum <- pat %>% group_by_all() %>% summarise(n = n())
    pat$which_pat <- match(pat$pat, pat_sum$pat)
    p <- ncol(S)
    z <- matrix(0, nrow = p, ncol = p)
    omega <- map(seq(nrow(pat_sum)), function(i){
      ixT <- which(pat_sum[i,] == TRUE)
      ixF <- which(pat_sum[i,] == FALSE)
      ii <- which(pat$which_pat == i)
      if(length(ixT) == 0){
        ome <- map(ii, function(j){
          s <- S[j,]
          solve_diag_psd_diag(R, s)
        })
        return(ome)
      }
      myR <- R[-ixT,-ixT, drop = FALSE]
      cat(dim(myR), "\n")
      pn <- nrow(myR)
      ome <- map(ii, function(j){
        s <- S[j,-ixT]
        o <- solve_diag_psd_diag(myR, diag(s, nrow = pn))
        om <- z
        om[ixF, ixF] <- o
        om
      })
      return(ome)
    }) %>% unlist(recursive = FALSE)
    ix1 <- sapply(seq(nrow(pat_sum)), function(i){which(pat$which_pat == i)}) %>% unlist()
    omega <- omega[match(seq(nrow(S)), ix1)]
  }
  return(omega)
}

get_omega_logdet <- function(omega, s_equal, n) {
  if(s_equal){
    # Log(det(omega))
    as.numeric(determinant(omega, logarithm = T)$modulus) * n
  }else{
    sum(
      sapply(omega, function(o) {
        as.numeric(determinant(o, logarithm = T)$modulus)
      })
    )
  }
}

set_data <- function(beta_hat_Y, se_Y, beta_hat_X, se_X, R,
                     ld_scores, RE, tau_init){

  beta_hat_X <- check_matrix(beta_hat_X)
  n <- nrow(beta_hat_X)
  p <- ncol(beta_hat_X)
  se_X <- check_matrix(se_X, n, p)
  if(!is.null(beta_hat_Y)){
    beta_hat_Y <- check_numeric(beta_hat_Y, n)
    se_Y <- check_numeric(se_Y, n)
    p <- p + 1
    beta_hat_X <- cbind(beta_hat_Y, beta_hat_X)
    se_X <- cbind(se_Y, se_X)
  }
  R <- check_matrix(R, p, p)
  R <- check_R(R)

  dat <- check_missing( beta_hat_X, se_X) # dat now has Y, S, s_equal, any_missing, n, and p
  dat$traits <- 1:p

  if(is.null(RE)){
    dat$omega <- get_omega(R, dat$S, dat$s_equal, dat$any_missing) # omega is row covariance of data, either list or single matrix
    # Pre-compute log(det(omega))
    dat$omega_logdet <- get_omega_logdet(dat$omega, dat$s_equal, n = dat$n)
    return(dat)
  }

  RE <- check_matrix(RE, p, p)
  dat$RE <- check_R(RE)
  dat$ld_scores <- check_numeric(ld_scores, n)

  dat$sigma <- get_sigma(R, dat$S, dat$s_equal, dat$any_missing)
  dat$tau <- tau_init
  dat$omega <- get_omega_tau(dat$sigma, dat$tau, dat$ld_scores, dat$RE)
  # Pre-compute log(det(omega))
  dat$omega_logdet <- get_omega_logdet(dat$omega, dat$s_equal, n = dat$n)
  dat$s_equal <- FALSE
  return(dat)
}


set_data_factors <- function(beta_hat_Y, se_Y, beta_hat_X, se_X,
                             beta_hat_Z, se_Z,
                             factors_matrix, factors_residual_sd,
                             R, ld_scores, RE, tau_init){

  if(is.null(beta_hat_Y)){
    stop("beta_hat_Y must be supplied for esmr_factors.\n")
  }
  beta_hat_Z <- check_matrix(beta_hat_Z)
  n <- nrow(beta_hat_Z)
  beta_hat_Y <- check_numeric(beta_hat_Y, n)
  beta_hat_X <- check_numeric(beta_hat_X, n)
  p <- ncol(beta_hat_Z) + 2

  se_Z <- check_matrix(se_Z, n, p-2)
  se_X <- check_numeric(se_X, n)
  se_Y <- check_numeric(se_Y, n)


  ## check factors
  factors_matrix <- check_matrix(factors_matrix, p-2 ) # F should be p-2 by k
  factors_residual_sd <- check_numeric(factors_residual_sd, p-2)
  k <- ncol(factors_matrix)

  if(!is.null(factors_residual_sd)){
    se_Z <- t(t(se_Z)*factors_residual_sd)
  }
  beta_hat_X <- cbind(beta_hat_X, beta_hat_Z)
  beta_hat_X <- cbind(beta_hat_Y, beta_hat_X)
  se_X <- cbind(se_X, se_Z)
  se_X <- cbind(se_Y, se_X)


  R <- check_matrix(R, p, p)
  R <- check_R(R)

  dat <- check_missing( beta_hat_X, se_X) # dat now has Y, S, s_equal, any_missing, n, and p
  dat$traits <- 1:p
  dat$factors_matrix <- factors_matrix
  dat$nfactors <- k
  dat$k <- k + 2


  if(is.null(RE)){
    dat$omega <- get_omega(R, dat$S, dat$s_equal, dat$any_missing) # omega is row covariance of data, either list or single matrix
    # Pre-compute log(det(omega))
    dat$omega_logdet <- get_omega_logdet(dat$omega, dat$s_equal, n = dat$n)
    return(dat)
  }

  RE <- check_matrix(RE, p, p)
  dat$RE <- check_R(RE)
  dat$ld_scores <- check_numeric(ld_scores, n)

  dat$sigma <- get_sigma(R, dat$S, dat$s_equal, dat$any_missing)
  dat$tau <- tau_init
  dat$omega <- get_omega_tau(dat$sigma, dat$tau, dat$ld_scores, dat$RE)
  # Pre-compute log(det(omega))
  dat$omega_logdet <- get_omega_logdet(dat$omega, dat$s_equal, n = dat$n)
  dat$s_equal <- FALSE


  return(dat)
}




order_upper_tri <- function(dat,
                            direct_effect_template,
                            direct_effect_init= NULL,
                            restrict_dag = TRUE){


  B <- direct_effect_template

  # Check if we have lower triangular
  if (any(B[upper.tri(B)] != 0) && restrict_dag) {
      # Direct effect template is not an lower triangular matrix
      # Attempt to re-order with topo-sort
      topo_order <- rlang::try_fetch({
        topo_sort_mat(B)
      }, error = function(cnd){
        rlang::abort(
          message = "Failed to find a lower triangular representation of the direct effect template. Check that supplied template corresponds to a valid DAG.\n",
          parent = cnd,
          call = rlang::call2("order_upper_tri"))
      })
      dat <- reorder_data(dat, topo_order)
      B <- B[topo_order, topo_order]
  }

  dat$B_template <- B
  if(!is.null(direct_effect_init)){
    o <- match(dat$traits, 1:dat$p)
    dat$B_init <- check_matrix(direct_effect_init, dat$p, dat$p)[o, o]
    if(any((dat$B_init != 0) & (dat$B_template == 0))) {
      rlang::abort("Initialization pattern does not match template.\n")
    }
  }else{
    dat$B_init <- matrix(0, nrow = dat$p, ncol = dat$p)
  }
  return(dat)
}



reorder_data <- function(
    dat, cols) {

  dat$Y <- dat$Y[,cols,drop=F]
  dat$S <- dat$S[,cols,drop=F]

  if(!is.null(dat$l)){
    dat$l$lbar <- dat$l$lbar[,cols,drop=F]
    dat$l$l2bar <- dat$l$l2bar[,cols,drop=F]
    #dat$l$abar <- dat$l$abar[,cols,drop=F]
    #dat$l$a2bar <- dat$l$a2bar[,cols,drop=F]
    dat$l$lfsr <- dat$l$lfsr[,cols,drop=F]
    dat$l$g_hat <- dat$l$g_hat[cols,drop=F]
  }

  if(!is.null(dat[["beta"]])) {
    dat$beta$beta_j <- match(dat$beta$beta_j, table = cols)
    dat$beta$beta_k <- match(dat$beta$beta_k, table = cols)
    dat$f <- make_f(dat)
  }
  if(!is.null(dat$omega)) {
    if(dat$s_equal){
      dat$omega <- dat$omega[cols, cols]
    }else{
      dat$omega <- lapply(dat$omega, function(x) x[cols, cols])
    }
  }
  if(!is.null(dat$G)){
    dat$G <- dat$G[cols,]
  }
  if(!is.null(dat$B_template)){
    dat$B_template <- dat$B_template[cols, cols]
  }
  if(!is.null(dat$B_init)){
    dat$B_init <- dat$B_init[cols,cols]
  }
  dat$traits <- dat$traits[cols]
  return(dat)
}

subset_data <- function(dat, ix){
  #s_equal <- check_equal_omega(dat$omega)
  dat$Y <- dat$Y[ix,,drop=F]
  dat$S <- dat$S[ix,,drop=F]
  dat$n <- length(ix)
  dat$l$lbar <- dat$l$lbar[ix,,drop=F]
  dat$l$l2bar <- dat$l$l2bar[ix,,drop=F]
  dat$l$abar <- dat$l$abar[ix,,drop=F]
  dat$l$a2bar <- dat$l$a2bar[ix,,drop=F]
  dat$l$lfsr <- dat$l$lfsr[ix,,drop=F]
  if(!dat$s_equal){
    dat$omega <- dat$omega[ix]
  }
  return(dat)
}

get_ix1_ix0 <- function(dat, ix1, remove_empty_B_cols = FALSE){
  if("integer" %in% class(ix1) | "numeric" %in% class(ix1)){
    stopifnot(all(ix1 %in% (1:dat$n)))
    dat$ix1 <- sort(ix1)
  }else if("character" %in% class(ix1)){
    ix1 <- stringr::str_split(ix1, "-", n = 2)[[1]]
    type <- ix1[1]
    thresh <- as.numeric(ix1[2])
    if (dat$is_nesmr && remove_empty_B_cols) {
      out_order <- rowSums(dat$B_template != 0)
      out_ix <- which(out_order > 0)
    } else if (dat$is_nesmr && !remove_empty_B_cols) {
      out_ix <- 1:dat$p
    } else {
      # Remove first column for esmr
      out_ix <- -1
    }

    if(type == "pval"){
      pval <- with(dat, 2*pnorm(-abs(Y/S)))
      # When out_ix is empty (e.g., template is all zeros), keep all variants
      if (length(out_ix) == 0) {
        dat$ix1 <- 1:dat$n
      } else {
        vals <- apply(pval[,out_ix,drop = FALSE], 1, min)
        dat$ix1 <- which(vals < thresh)
      }
    }else{
      stop("Unknown option to ix1\n")
    }
  }else{
    stop("Unknown option to ix1\n")
  }
  dat$ix0 <- setdiff((1:dat$n), dat$ix1) |> sort()
  return(dat)
}

delta_method_pvals <- function(dat){
  e_ix <- which(!dat$beta$fix_beta)
  fix_ix <- which(dat$beta$fix_beta)
  e_coords <- cbind(dat$beta$beta_k, dat$beta$beta_j)[e_ix,,drop=FALSE]
  if(length(fix_ix) > 0){
    fix_coords <- cbind(dat$beta$beta_k, dat$beta$beta_j)[fix_ix,,drop=FALSE]
    colnames(fix_coords) <- c("row", "col")
  }

  f <- function(tot){
    myT <- matrix(0, nrow = dat$p, ncol = dat$p)
    myT[e_coords] <- tot

    if(length(fix_ix) > 0){
      myT <- complete_T(myT, fix_coords)$total_effects
    }

    myB <- total_to_direct(myT)
    dir <- myB[e_coords]
    return(dir)
  }
  jac <- numDeriv::jacobian(f, x = dat$beta$beta_m[e_ix])
  V <- dat$beta$V[e_ix,e_ix]
  VB <- jac %*% V %*% t(jac)
  muB <- f(dat$beta$beta_m[e_ix])
  log_pvals <- log(2) + pnorm(-abs(muB/sqrt(diag(VB))), log.p = TRUE)
  pmat <- semat <- matrix(0, nrow = dat$p, ncol = dat$p)
  pmat[e_coords] <- log_pvals
  semat[e_coords] <- sqrt(diag(VB))
  return(list(pmat = pmat, semat = semat))
}

## unused
get_wpost <- function(beta_hat, se_beta_hat, col_ix, prior_family = "point_normal"){
  wpost <- purrr::map_dfc(col_ix, function(ii){
    cat(ii, "\n")
    x <- beta_hat[,ii];
    s <- se_beta_hat[,ii];
    f <- ebnm(x = x, s = s, prior_family = "point_normal", output = ebnm::output_all());
    pi0 <- f$fitted_g$pi[1];
    mu <- f$fitted_g$mean[2];
    s2 <- f$fitted_g$sd[2]^2;
    w <- 1-pi0;
    a <- 1/s2;
    wpost <- ebnm:::wpost_normal(x=x, s=s, w, a, mu);
    df <- data.frame(w = wpost);
    names(df) <- paste0("wpost", ii);
    return(df);
  });
  return(wpost)
}

get_lower_triangular <- function(x, diag = FALSE) {
  x[lower.tri(x, diag)]
}

# Converts a matrix to an edgelist: from->to with the value as the matrix value
matrix_to_edgelist <- function(
    X, lower_tri = FALSE, value = 'value',
    remove_diag = FALSE) {
  if (lower_tri) {
    ltx <- lower.tri(X, diag = ! remove_diag)
  } else {
    ltx <- TRUE
  }
  res <- data.frame(
    from = row(X)[ltx],
    to = col(X)[ltx]
  )
  res[[value]] <- X[ltx]
  if (remove_diag) {
    res <- res[res$from != res$to, ]
  }
  res
}

flat_string_to_adj_mat <- function(x) {
  # First get a vector of the digits of x
  digits <- as.numeric(unlist(strsplit(as.character(x), "")))
  # Then turn into a matrix
  n <- sqrt(length(digits))
  stopifnot(n == floor(n))
  matrix(digits, nrow = n, ncol = n)
}
