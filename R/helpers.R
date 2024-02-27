default_precision <- function(dims){
  sqrt(.Machine$double.eps)*prod(dims)
}

make_f <- function(dat){
  if(length(dat$beta$beta_j) == 0 | is.null(dat$beta$beta_j)){
    return(list(fbar = diag(dat$p), f2bar = diag(dat$p),
                fgbar = dat$G, fg2bar = dat$G^2))
  }
  fbar <- f2bar <- diag(dat$p)
  nb <- length(dat$beta$beta_j)
  ix <- cbind(dat$beta$beta_j, dat$beta$beta_k)
  fbar[ix] <- dat$beta$beta_m

  fgbar <- fbar %*% dat$G
  V <- matrix(0, nrow = dat$p, ncol = dat$p)
  V[ix] <- dat$beta$beta_s^2

  f2bar <- (fbar^2) + V
  fg2bar <- (fgbar^2) + (V %*% (dat$G^2))
  return(list(fgbar = fgbar, fg2bar = fg2bar,
              fbar = fbar, f2bar = f2bar))
}

get_omega <- function(R, S, any_missing){

  s_equal <- apply(S, 2, function(x){all(x == x[1])}) %>% all()
  p <- ncol(S)
  R_is_id <- is.null(R) | all(R == diag(p))

  if(s_equal & R_is_id){
    s <- S[1,]
    omega <- diag(1/s^2)
  }else if(s_equal){
    s <- S[1,]
    omega <- solve_diag_psd_diag(R, s)
  }else if(R_is_id & !any_missing){
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

set_data <- function(beta_hat_Y, se_Y, beta_hat_X, se_X, R){

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

  dat <- check_missing( beta_hat_X, se_X)
  dat$omega <- get_omega(R, dat$S, dat$any_missing) # omega is row correlation of data, either list or single matrix
  dat$n <- n
  dat$p <- p
  dat$traits <- 1:p
  return(dat)

}

order_upper_tri <- function(dat, direct_effect_template = NULL, direct_effect_init= NULL){
  if(!is.null(direct_effect_template)){
    B <- direct_effect_template
  }else{
    B <- matrix(0, nrow = dat$p, ncol = dat$p)
    B[2:dat$p, 1] <- 1
  }
  # Check if we have lower triangular
  if (any(B[upper.tri(B)] != 0)) {
      # Direct effect template is not an lower triangular matrix
      # Attempt to re-order with topo-sort
      topo_order <- tryCatch({
        topo_sort_mat(B)
      }, error = function(e){
        stop("Failed to find a lower triangular representation of the direct effect template. Check that supplied template corresponds to a valid DAG.\n")
      })
      dat <- reorder_data(dat, topo_order)
      # beta_hat_X <- beta_hat_X[, topo_order]
      # se_X <- se_X[, topo_order]
      # R <- R[topo_order, topo_order]
      B <- B[topo_order, topo_order]
  }

  dat$B_template <- B
  if(!is.null(direct_effect_init)){
    o <- match(dat$traits, 1:dat$p)
    dat$B_init <- check_matrix(direct_effect_init, dat$p, dat$p)[o, o]
    if(any(dat$B_init[!dat$B_template == 0] != 0)){
      stop("Initialization pattern does not match template.\n")
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
    dat$l$abar <- dat$l$abar[,cols,drop=F]
    dat$l$a2bar <- dat$l$a2bar[,cols,drop=F]
    dat$l$lfsr <- dat$l$lfsr[,cols,drop=F]
    dat$l$g_hat <- dat$l$g_hat[cols,drop=F]
  }
  if(!is.null(dat$f)){
    dat$f <- lapply(dat$f, function(x) {
      x[cols, cols]
    })
  }
  if(!is.null(dat$beta)){
    dat$beta$beta_j <- match(dat$beta$beta_j, table = cols)
    dat$beta$beta_k <- match(dat$beta$beta_k, table = cols)
  }
  if(!is.null(dat$omega)) {
    dat$omega <- lapply(dat$omega, function(x) x[cols, cols])
  }
  if(!is.null(dat$G)){
    dat$G <- dat$G[cols,cols]
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
  s_equal <- check_equal_omega(dat$omega)
  dat$Y <- dat$Y[ix,,drop=F]
  dat$S <- dat$S[ix,,drop=F]
  dat$n <- length(ix)
  dat$l$lbar <- dat$l$lbar[ix,,drop=F]
  dat$l$l2bar <- dat$l$l2bar[ix,,drop=F]
  dat$l$abar <- dat$l$abar[ix,,drop=F]
  dat$l$a2bar <- dat$l$a2bar[ix,,drop=F]
  dat$l$lfsr <- dat$l$lfsr[ix,,drop=F]
  if(!s_equal){
    dat$omega <- dat$omega[ix]
  }
  return(dat)
}

get_ix1_ix0 <- function(dat, ix1){
  if("integer" %in% class(ix1) | "numeric" %in% class(ix1)){
    stopifnot(all(ix1 %in% (1:dat$n)))
    dat$ix1 <- sort(ix1)
  }else if("character" %in% class(ix1)){
    ix1 <- stringr::str_split(ix1, "-", n = 2)[[1]]
    type <- ix1[1]
    thresh <- as.numeric(ix1[2])
    out_order <- rowSums(dat$B_template != 0)
    out_ix <- which(out_order > 0)
    if(type == "pval"){
      pval <- with(dat, 2*pnorm(-abs(Y/S)))
      vals <- apply(pval[,out_ix,drop = FALSE], 1, min)
      dat$ix1 <- which(vals < thresh)
    }else{
      stop("Unknown option to ix1\n")
    }
  }else{
    stop("Unknown option to ix1\n")
  }
  dat$ix0 <- setdiff((1:dat$n), dat$ix1) |> sort()
  return(dat)
}


## check structure of direct effects, return structure of total effects
direct_to_total <- function(B_dir){
  n <- nrow(B_dir)
  B_total <- solve(diag(n) - B_dir) - diag(n)
  if(!all(diag(B_total) == 0)){
    stop("Failed to compute total effects from direct. Check that supplied B_dir corresponds to a valid DAG.\n")
  }
  return(B_total)
}

total_to_direct <- function(B_tot){
  n <- nrow(B_tot)
  B_dir <-diag(n) - solve(diag(n) + B_tot)
  if(!all(diag(B_dir) == 0)){
    stop("Failed to compute total effects from direct. Check that supplied B_tot corresponds to a valid DAG.\n")
  }
  return(B_dir)
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
