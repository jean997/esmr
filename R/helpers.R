default_precision <- function(dims){
  sqrt(.Machine$double.eps)*prod(dims)
}

make_fbar <- function(beta_m, beta_s, beta_j, beta_k, p){
  fbar <- f2bar <- diag(p)
  for(i in seq(length(beta_j))) fbar[beta_j[i],beta_k[i]] <- beta_m[i]
  for(i in seq(length(beta_j))) f2bar[beta_j[i],beta_k[i]] <- beta_m[i]^2 + beta_s[i]^2
  return(list(fbar = fbar, f2bar = f2bar))
}

make_f_fun <- function(p, beta_j, beta_k){
  if(length(beta_j) == 0 | is.null(beta_j)){
    ff <- function(beta_m, beta_s){
      return(list(fbar = diag(p), f2bar = diag(p)))
    }
    return(ff)
  }
  ff <- function(beta_m, beta_s){
    fbar <- f2bar <- diag(p)
    for(i in seq(length(beta_j))) fbar[beta_j[i],beta_k[i]] <- beta_m[i]
    for(i in seq(length(beta_j))) f2bar[beta_j[i],beta_k[i]] <- beta_m[i]^2 + beta_s[i]^2
    return(list(fbar = fbar, f2bar = f2bar))
  }
  return(ff)
}

init_beta <- function(p, which_beta=NULL, beta_joint = TRUE, beta_m_init = NULL){
  if(!is.null(which_beta) & beta_joint){
    stop("Joint update only implemented for standard problem.\n")
  }
  if(is.null(which_beta) | beta_joint){
    beta_j <- rep(1, p-1)
    beta_k <- 2:p
  }else if(!is.null(which_beta)){
    beta_j <- which_beta[,1]
    beta_k <- which_beta[,2]
  }
  el <- length(beta_j)
  if(is.null(beta_m_init)){
    beta_m <- rep(0, el)
    beta_s <- rep(0, el)
  }else{
    if(!length(beta_m_init) == el) stop(paste0("Expected beta_m_init to have length ", el, ", found ", length(beta_m_init), "\n"))
    beta_m <- beta_m_init
    beta_s <- rep(0, el)
  }
  return(list(beta_m = beta_m, beta_s = beta_s, beta_j = beta_j, beta_k = beta_k))
}

init_l <- function(n, p){
  lbar <- l2bar <- matrix(0, nrow = n, ncol = p)
  lfsr <- matrix(1, nrow = n, ncol = p)
  g_hat <- list()
  return(list(lbar = lbar, l2bar = l2bar, lfsr = lfsr, g_hat = g_hat))
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
  R <- check_matrix(R, "se_X", p, p)
  R <- check_R(R)

  dat <- check_missing(cbind(beta_hat_Y, beta_hat_X), cbind(se_Y, se_X))
  dat$omega <- get_omega(R, dat$S, dat$any_missing)
  dat$n <- n
  dat$p <- p

  return(dat)

}
