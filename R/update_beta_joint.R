
#'@export
update_beta_joint <- function(Y, lbar, l2bar, omega, prior_cov = NULL){
  n <- nrow(Y)
  p <- ncol(Y)
  stopifnot(nrow(lbar) == n & ncol(lbar) == p)
  stopifnot(nrow(l2bar) == n & ncol(l2bar) == p)

  s_equal <- check_equal_omega(omega)

  if(is.null(prior_cov)){
    T0 <- matrix(0, nrow = p-1, ncol = p-1)
  }else{
    T0 <- check_matrix(prior_cov, p-1, p-1)
    T0 <- solve(prior_cov)
  }
  s2l <- l2bar - (lbar^2)
  if(s_equal){
    ltl <- t(lbar) %*% lbar
    diag(ltl) <- colSums(l2bar)
    R <- omega[1,1]*(ltl[-1,-1, drop = FALSE])
    R <- R + T0
    G <- t(t(Y - lbar)*omega[1,])
    g <- rowSums(G)
    a <- colSums(g*lbar[,-1,drop=FALSE]) - rowSums(t(s2l[,-1,drop=FALSE])*omega[1,-1])
    a <-  matrix(a, nrow = p-1)
    S <- solve(R)
    mu <- S %*% a
  }else{
    o11 <- map(omega, function(o){o[1,1]}) %>% unlist()
    O11 <- map(omega, function(o){o[1,]}) %>% unlist() %>%
      matrix(nrow = n, byrow = TRUE)
    ltl <- t(lbar) %*% (lbar*o11)
    diag(ltl) <- colSums(l2bar*o11)
    R <- ltl[-1,-1, drop = FALSE] + T0
    G <- (Y - lbar)*O11
    g <- rowSums(G)
    a <- colSums(g*lbar[,-1,drop = FALSE]) - colSums(s2l[,-1,drop = FALSE]*O11[,-1,drop=FALSE])
    a <-  matrix(a, nrow = p-1)
    S <- solve(R)
    mu <- S %*% a

  }
  return(list(m = mu, S = S))
}
