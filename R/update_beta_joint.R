
#'@export
update_beta_joint <- function(dat, j=1, ix = NULL, prior_cov = NULL, return_W = FALSE){

  p <- dat$p
  n <- dat$n
  if(is.null(ix)){
    ix <- seq(p)[-j]
  }else{
    stopifnot(all(ix %in% seq(p)))
    stopifnot(!any(duplicated(ix)))
    #ix <- sort(ix)
  }
  m <- length(ix)
  if(is.null(prior_cov)){
    T0 <- matrix(0, nrow = m, ncol = m)
  }else{
    T0 <- check_matrix(prior_cov, m, m)
    T0 <- solve(prior_cov)
  }
  Va <- dat$l$a2bar - (dat$l$abar^2)

  if(dat$s_equal){
    A <- t(dat$l$abar) %*% dat$l$abar + diag(colSums(Va))
    Astar <- dat$G %*% A %*% t(dat$G)

    Rfull <- dat$omega[j,j]*Astar  # W in the manuscript
    a10 <- colSums(dat$l$lbar *rowSums(t(t(dat$Y)*dat$omega[,j])))
    a20 <- lapply(seq(p)[-j], function(jj){
      Astar%*% t(dat$f$fbar[jj,,drop = FALSE])*dat$omega[j,jj]
    }) %>% Reduce(`+`, .)
    afull <- matrix(a10 - a20, nrow = p)
  }else{
    Oj <- map(dat$omega, function(o){o[j,]}) %>% unlist() %>%
      matrix(nrow = n, byrow = TRUE)
    Astar <- lapply(seq(p), function(jj){ # this is a list of W^{(a,j)}
      A <- t(dat$l$abar * Oj[,jj]) %*% dat$l$abar + diag(colSums(Va * Oj[,jj]))
      dat$G %*% A %*% t(dat$G)
    })
    Rfull <- Astar[[j]]
    a10 <- colSums(dat$l$lbar *rowSums(dat$Y*Oj))
    a20 <- lapply(seq(p)[-j], function(jj){
      Astar[[jj]]%*% t(dat$f$fbar[jj,,drop = FALSE])
    }) %>% Reduce(`+`, .)
    afull <- matrix(a10 - a20, nrow = p)
  }
  if(length(ix) < p){
    R <- Rfull[ix,ix]
    R12 <- Rfull[ix,-ix]
    a <- afull[ix] - R12 %*% t(dat$f$fbar[j,-ix,drop = FALSE])
  }else{
    R <- Rfull
    a <- afull
  }

  tryCatch({
    S <- solve(R)
  }, error = function(e){
    all_zero_cols_lbar <- apply(zapsmall(dat$l$lbar, digits = 10), 2, function(x){
      all(x == 0)
    })
    if (any(all_zero_cols_lbar)) {
      stop('lbar has a column of all zeros for column(s): ', which(all_zero_cols_lbar))
    } else {
      stop('Cannot invert R matrix: ', R)
    }
  })
  mu <- S %*% a

  if(return_W){
    return(list(m = mu, S = S, W = R, b = a))
  }
  return(list(m = mu, S = S))
}


update_beta_full_joint <- function(dat, prior_cov = NULL){

  p <- dat$p
  n <- dat$n
  k <- dat$k

  ix <- p*(dat$beta$beta_k-1) + dat$beta$beta_j
  ix <- ix[!dat$beta$fix_beta]
  m <- length(ix)

  if(is.null(prior_cov)){
    T0 <- matrix(0, nrow = m, ncol = m)
  }else{
    T0 <- check_matrix(prior_cov, m, m)
    T0 <- solve(prior_cov)
  }

  Va <- dat$l$a2bar - (dat$l$abar^2)

  if(dat$s_equal){
    A <- t(dat$l$abar) %*% dat$l$abar + diag(colSums(Va))
    Astar <- dat$G %*% A %*% t(dat$G)
    Rfull <- kronecker(Astar, dat$omega)

    OYt <- dat$omega %*% t(dat$Y)
    afull <- lapply(seq(n), function(i){
      kronecker( matrix(dat$l$lbar[i,], nrow = k), matrix(OYt[,i], nrow = p))
    }) %>% Reduce(`+`, .)
  }else{
    Rfull <- lapply(seq(n), function(i){
      l <- dat$l$abar[i,]
      a <- outer(l, l) + diag(Va[i,], nrow = k)
      kronecker(tcrossprod(dat$G, tcrossprod(dat$G, a)), dat$omega[[i]]) ## kronecker(G %*% a %*% t(G), O)
    }) %>% Reduce(`+`, .)
    afull <- lapply(seq(n), function(i){
      kronecker(dat$l$lbar[i,], dat$omega[[i]] %*% matrix(dat$Y[i,], nrow = p))
    }) %>% Reduce(`+`, .)
  }
  if(length(ix) < p*k){
    R <- Rfull[ix,ix]
    R12 <- Rfull[ix,-ix]
    fb <- matrix(as.vector(dat$f$fbar), ncol = 1)
    a <- afull[ix] - R12 %*% fb[-ix,,drop = F]
  }else{
    R <- Rfull
    a <- afull
  }
  S <- solve(R)
  mu <- S %*% a
  return(list(m = mu, S = S))
}
