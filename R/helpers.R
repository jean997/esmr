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
  for(i in seq(nb)) fbar[dat$beta$beta_j[i],dat$beta$beta_k[i]] <- dat$beta$beta_m[i]
  fgbar <- fbar %*% dat$G
  V <- matrix(0, nrow = dat$p, ncol = dat$p)
  for(i in seq(nb)) V[dat$beta$beta_j[i],dat$beta$beta_k[i]] <- dat$beta$beta_s[i]^2
  f2bar <- (fbar^2) + V
  fg2bar <- (fgbar^2) + (V %*% (dat$G^2))
  return(list(fgbar = fgbar, fg2bar = fg2bar,
              fbar = fbar, f2bar = f2bar))
}


init_beta <- function(p, which_beta=NULL,
                      beta_m_init = NULL,
                      fix_beta = FALSE){
  #if(!is.null(which_beta) & beta_joint){
  #  stop("Joint update only implemented for standard problem.\n")
  #}
  if(is.null(which_beta)){
    beta_j <- rep(1, p-1)
    beta_k <- 2:p
  }else if(!is.null(which_beta)){
    beta_j <- which_beta[,1]
    beta_k <- which_beta[,2]
  }
  nb <- length(beta_j)
  if(is.null(beta_m_init)){
    beta_m <- rep(0, nb)
    beta_s <- rep(0, nb)
  }else{
    if(!length(beta_m_init) == nb) stop(paste0("Expected beta_m_init to have length ", nb, ", found ", length(beta_m_init), "\n"))
    beta_m <- beta_m_init
    beta_s <- rep(0, nb)
  }

  if(length(fix_beta) == 1){
    fix_beta <- rep(fix_beta, nb)
  }else if(length(fix_beta) != nb){
    stop("fix_beta and expected to have length 1 or ", nb, ". Found ", length(fix_beta), ".\n")
  }
  stopifnot(class(fix_beta) == "logical")
  return(list(beta_m = beta_m, beta_s = beta_s,
              beta_j = beta_j, beta_k = beta_k, fix_beta = fix_beta))
}


# previously init_l_future
init_l <- function(n, p, m){
  lbar <- matrix(0, nrow = n, ncol = p)
  abar <- matrix(0, nrow = n, ncol = m)
  lfsr <- matrix(1, nrow = n, ncol = m)
  g_hat <- list()
  return(list(lbar = lbar, l2bar = lbar,
              abar = abar, a2bar = abar,
              #wpost = lbar, mupost = lbar, s2post = lbar,
              #post_mode = abar,
              lfsr = lfsr, g_hat = g_hat))
}


check_numeric <- function(x, string, n){
  if(is.null(x)) return(x)
  if("matrix" %in% class(x) | "data.frame" %in% class(x)){
    if(ncol(x) > 1) stop(paste0(string, " must be a numeric vector or one column array."))
    x <- as.numeric(x[,1])
  }else if(!"numeric" %in% class(x)){
    stop(paste0(string, " must be a numeric vector or one column array."))
  }
  if(!missing(n)){
     if(length(x) != n) stop(paste0("Expected ", string, " to have length ", n, ", found ", length(x), "\n"))
  }
  return(x)
}

check_matrix <- function(x, string, n, p){
  if(is.null(x)) return(x)
  if("data.frame" %in% class(x)){
    cat("Coercing ", string, " to matrix.\n")
    x <- as.matrix(x)
  }else if("numeric" %in% class(x)){
    x <- as.matrix(x, ncol = 1)
  }else if(!"matrix" %in% class(x)){
    stop(paste0(string, " must be a numeric vector, matrix, or data.frame."))
  }
  if(!missing(n)){
    if(nrow(x) != n) stop(paste0("Expected ", string, " to have ", n, " rows, found ", nrow(x), "\n"))
  }
  if(!missing(p)){
    if(ncol(x) != p) stop(paste0("Expected ", string, " to have ", p, " columns, found ", ncol(x), "\n"))
  }
  return(x)
}

check_R <- function(R, tol = 1e-8){
  if(is.null(R)) return(R)
  if(!Matrix::isSymmetric(R)){
    stop("R is not symmetric.\n")
  }
  evR <- eigen(R, only.values = TRUE)$values
  if(!all(evR > tol)){
    stop("R is not positive definite.\n")
  }
  if(!all(diag(R) == 1)){
    stop("R should be a correlation matrix.\n")
  }
  return(R)
}

check_missing <- function(Y, S){
  missing_ix <- which(is.na(Y))
  any_missing <- length(missing_ix) > 0
  if(length(missing_ix) == 0){
    return(list(Y = Y, S = S, any_missing = FALSE))
  }
  if(any(is.na(S[-missing_ix]))) stop("Found missing SEs for non-missing effect estimates.\n")

  nmiss_r <- rowSums(is.na(Y))
  nmiss_c <- colSums(is.na(Y))
  p <- ncol(Y)
  n <- nrow(Y)
  if(any(nmiss_r == p)) stop("Data cannot have variants missing for all traits.\n")
  if(any(nmiss_c == n)) stop("At least one trait is missing for all variants.\n")
  Y[missing_ix] <- 0
  S[missing_ix] <- NA
  return(list(Y = Y, S = S, any_missing = TRUE))
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
    omega <- solve(diag(s) %*% R %*% diag(s))
  }else if(R_is_id & !any_missing){
    omega <- apply(S, 1, function(s){
      diag(1/s^2, nrow = p)
    }, simplify = FALSE)
  }else if(!any_missing){
    omega <- apply(S, 1, function(s){
      solve(diag(s) %*% R %*% diag(s))
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
          solve(diag(s) %*% R %*% diag(s))
        })
        return(ome)
      }
      myR <- R[-ixT,-ixT, drop = FALSE]
      cat(dim(myR), "\n")
      pn <- nrow(myR)
      ome <- map(ii, function(j){
        s <- S[j,-ixT]
        o <- solve(diag(s, nrow = pn) %*% myR %*% diag(s, nrow = pn))
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

  beta_hat_Y <- check_numeric(beta_hat_Y, "beta_hat_Y")
  n <- length(beta_hat_Y)
  se_Y <- check_numeric(se_Y, "se_Y", n)
  beta_hat_X <- check_matrix(beta_hat_X, "beta_hat_X", n)
  p <- ncol(beta_hat_X) + 1
  se_X <- check_matrix(se_X, "se_X", n, p-1)
  R <- check_matrix(R, "R", p, p)
  R <- check_R(R)

  dat <- check_missing(cbind(beta_hat_Y, beta_hat_X), cbind(se_Y, se_X))
  dat$omega <- get_omega(R, dat$S, dat$any_missing) # omega is row correlation of data, either list or single matrix
  dat$n <- n
  dat$p <- p

  return(dat)

}

subset_data <- function(dat, ix){
  s_equal <- check_equal_omega(dat$omega)
  dat$Y <- dat$Y[ix,]
  dat$S <- dat$S[ix,]
  dat$n <- length(ix)
  dat$l$lbar <- dat$l$lbar[ix,]
  dat$l$l2bar <- dat$l$l2bar[ix,]
  dat$l$abar <- dat$l$abar[ix,]
  dat$l$a2bar <- dat$l$a2bar[ix,]
  dat$l$lfsr <- dat$l$lfsr[ix,]
  if(!s_equal){
    dat$omega <- dat$omega[ix]
  }
  return(dat)
}

check_equal_omega <- function(omega){
  if("matrix" %in% class(omega)){
    #check_matrix(omega, "omega", p, p)
    s_equal <- TRUE
  }else{
    stopifnot(class(omega) == "list")
    s_equal <- FALSE
  }
  return(s_equal)
}

get_ix1_ix0 <- function(dat, ix1){
  if("integer" %in% class(ix1) | "numeric" %in% class(ix1)){
    stopifnot(all(ix1 %in% (1:dat$n)))
    dat$ix1 <- sort(ix1)
  }else if("character" %in% class(ix1)){
    ix1 <- stringr::str_split(ix1, "-", n = 2)[[1]]
    type <- ix1[1]
    thresh <- as.numeric(ix1[2])
    if(type == "pval"){
      pval <- with(dat, 2*pnorm(-abs(Y/S)))
      vals <- apply(pval[,-1,drop = FALSE], 1, min)
      dat$ix1 <- which(vals < thresh)
    }else if(type == "lfsr" | type == "zl"){
      for(i in 1:5){
        dat <- update_l_sequential(dat)
      }
      if(type == "lfsr"){
        vals <- apply(dat$l$lfsr[,-1, drop = F], 1, min)
        dat$ix1 <- which(vals < thresh)
      }else{
        vl <- with(dat$l, l2bar - (lbar^2))
        zl <- dat$l$lbar/sqrt(vl)
        vals <- apply(abs(zl[,-1,drop = F]), 1,max)
        dat$ix1 <- which(vals > thresh)
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

## unused
check_omega <- function(omega, n, p, s_equal){
  if(s_equal){
    if(!"matrix" %in% class(omega)) stop("omega is not of class matrix\n")
    check_matrix(omega, "omega", p, p)
  }else{
    stopifnot(class(omega) == "list")
    stopifnot(length(omega) == n)
    d <- lapply(omega, dim)
    nr <- map(d, 1) %>% unlist()
    nc <- map(d, 2) %>% unlist()
    stopifnot(all(nr == p))
    stopifnot(all(nc == p))
  }
}


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
