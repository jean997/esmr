
#'@export
update_beta_joint <- function(dat, j=1, ix = NULL, prior_cov = NULL){
  #stopifnot(j == 1)
  #Y, lbar, l2bar, omega, prior_cov = NULL){
  #n <- nrow(Y)
  #p <- ncol(Y)
  #stopifnot(nrow(lbar) == n & ncol(lbar) == p)
  #stopifnot(nrow(l2bar) == n & ncol(l2bar) == p)

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
  #s2l <- l2bar - (lbar^2)
  Vl <- dat$l$l2bar_o - (dat$l$lbar_o^2)

  if(s_equal){
    A <- t(dat$l$lbar_o) %*% dat$l$lbar_o + diag(colSums(Vl))
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
      Astar%*% dat$f$fbar_o[,jj,drop = FALSE]*dat$omega[j,jj]
    }) %>% Reduce(`+`, .)
    afull <- matrix(a10 - a20, nrow = p)
  }else{
    Oj <- map(dat$omega, function(o){o[j,]}) %>% unlist() %>%
      matrix(nrow = n, byrow = TRUE)
    Astar <- lapply(seq(p), function(jj){
      A <- t(dat$l$lbar_o * Oj[,jj]) %*% dat$l$lbar_o + diag(colSums(Vl * Oj[,jj]))
      dat$G %*% A %*% t(dat$G)
    })
    Rfull <- Astar[[j]]
    a10 <- colSums(dat$l$lbar *rowSums(dat$Y*Oj))
    a20 <- lapply(seq(p)[-j], function(jj){
      Astar[[jj]]%*% dat$f$fbar_o[,jj,drop = FALSE]
    }) %>% Reduce(`+`, .)
    afull <- matrix(a10 - a20, nrow = p)
  }
  if(length(ix) < p){
    R <- Rfull[ix,ix]
    R12 <- Rfull[ix,-ix]
    a <- afull[ix] - R12 %*% t(dat$f$fbar_o[j,-ix,drop = FALSE])
  }else{
    R <- Rfull
    a <- afull
  }
  S <- solve(R)
  mu <- S %*% a
  return(list(m = mu, S = S))
}
