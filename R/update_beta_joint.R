
#'@export
update_beta_joint <- function(dat, j=1, ix = NULL, prior_cov = NULL, return_W = FALSE){
  s_equal <- check_equal_omega(dat$omega)

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

  if(s_equal){
    A <- t(dat$l$abar) %*% dat$l$abar + diag(colSums(Va))
    Astar <- dat$G %*% A %*% t(dat$G)
    # R <- dat$omega[j,j]*Astar[ix,ix] + T0
    # a1 <- colSums(dat$l$lbar[,ix]*rowSums(t(t(dat$Y)*dat$omega[,j])))
    # a2 <- colSums(Astar[,ix]*dat$omega[,j])
    # a <- matrix(a1 - a2, nrow = m)
    # S <- solve(R)
    # mu <- S %*% a

    Rfull <- dat$omega[j,j]*Astar
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
  S <- solve(R)
  mu <- S %*% a

  if(return_W){
    return(list(m = mu, S = S, W = R, b = a))
  }
  return(list(m = mu, S = S))
}


update_beta_full_joint <- function(dat, prior_cov = NULL){
  s_equal <- check_equal_omega(dat$omega)

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

  if(s_equal){
    A <- t(dat$l$abar) %*% dat$l$abar + diag(colSums(Va))
    Astar <- dat$G %*% A %*% t(dat$G)
    Rfull <- kronecker(Astar, dat$omega)

    OYt <- dat$omega %*% t(dat$Y)
    afull <- lapply(seq(n), function(i){
      kronecker( matrix(dat$l$lbar[i,], nrow = k), matrix(OYt[,i], nrow = p))
    }) %>% Reduce(`+`, .)
  }else{
    Rfull <- lapply(seq(n), function(i){
      l <- matrix(dat$l$abar[i,], nrow = k)
      a <- l %*% t(l) + diag(Va[i,], nrow = k)
      kronecker(tcrossprod(dat$G, tcrossprod(dat$G, a)), dat$omega[[i]]) ## kronecker(G %*% a %*% t(G), O)
    }) %>% Reduce(`+`, .)
    afull <- lapply(seq(n), function(i){
      kronecker(matrix(dat$l$lbar[i,], nrow = k), dat$omega[[i]] %*% matrix(dat$Y[i,], nrow = p))
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
