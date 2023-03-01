
# Calculates E_q[log lik(l, f | Y, omega)]
# The objective function is E_q[log lik(l, f | Y, omega)] - KL(q_l || g_l)
# The second term (including minus sign) is returned as res$obj (confusingly named)
# so elbo = calc_ell + res$obj
# need to make this calculation faster
# Should maybe take f2bar out and treat fbar as fixed
#'@export
#calc_ell <- function(Y, lbar, l2bar, fbar, f2bar, omega){
calc_ell <- function(Y, lbar, l2bar, fbar, omega){
  n <- nrow(Y)
  p <- ncol(Y)
  stopifnot(nrow(lbar) == n & ncol(lbar) == p)
  stopifnot(nrow(l2bar) == n & ncol(l2bar) == p)
  #stopifnot(nrow(fbar) == p & nrow(f2bar) == p & ncol(fbar) == p & ncol(f2bar) == p)
  stopifnot(nrow(fbar) == p  & ncol(fbar) == p )


  if("matrix" %in% class(omega)){
    s_equal <- TRUE
  }else{
    stopifnot(class(omega) == "list")
    stopifnot(length(omega) == n)
    s_equal <- FALSE
  }

  ybar <- lbar %*% t(fbar)

  varlbar <- l2bar - (lbar^2)
  #varfbar <- f2bar - (fbar^2)
  if(s_equal){
    part_a <- map_dbl(seq(nrow(ybar)), function(i){
      crossprod(t(Y[i,,drop =F]), omega) %>% tcrossprod(ybar[i, , drop = F])
    }) %>% sum()
    part_b <- map_dbl(seq(nrow(ybar)), function(i){
      crossprod(t(ybar[i,,drop =F]), omega) %>% tcrossprod(ybar[i, , drop = F])
    }) %>% sum()
    #x1 <- sum(t( t(tcrossprod(lbar^2 , varfbar))*diag(omega)))
    x2 <- sum(t( t(tcrossprod(varlbar , fbar^2))*diag(omega)))
    #x3 <- sum(t( t(tcrossprod(varlbar , varfbar))*diag(omega)))
    #part_b <- part_b + x1 + x2 + x3
    part_b <- part_b + x2
    ell <- part_a - 0.5*part_b
  }else{
    part_a <- map_dbl(seq(nrow(ybar)), function(i){
      crossprod(t(Y[i,,drop =F]), omega[[i]]) %>% tcrossprod(ybar[i, , drop = F])
    }) %>% sum()
    part_b <- map_dbl(seq(nrow(ybar)), function(i){
      crossprod(t(ybar[i,,drop =F]), omega[[i]]) %>% tcrossprod(ybar[i, , drop = F])
    }) %>% sum()
    #x1 <- map_dbl(seq(nrow(ybar)), function(i){
    #  sum(t(t( tcrossprod((lbar^2)[i,,drop =F], varfbar))*diag(omega[[i]])) )
    #}) %>% sum()
    x2 <- map_dbl(seq(nrow(ybar)), function(i){
      sum(t(t( tcrossprod(varlbar[i,,drop =F], fbar^2))*diag(omega[[i]])) )
    }) %>% sum()
    #x3 <- map_dbl(seq(nrow(ybar)), function(i){
    #  sum(t(t( tcrossprod(varlbar[i,,drop =F], varfbar))*diag(omega[[i]])) )
    #}) %>% sum()
    #part_b <- part_b + x1 + x2 + x3
    part_b <- part_b +  x2
    ell <- part_a - 0.5*part_b
  }
  return(ell)
}


#'@export
calc_ell2 <- function(Y, lbar, l2bar, fbar, omega){
  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(fbar)
  check_matrix(lbar, "lbar", n, k)
  check_matrix(l2bar, "l2bar", n, k)
  check_matrix(fbar, "fbar", p, k)

  s_equal <- check_equal_omega(omega)

  ybar <- lbar %*% t(fbar)
  R <- Y - ybar

  varlbar <- l2bar - (lbar^2)
  #varfbar <- f2bar - (fbar^2)
  if(s_equal){
    part_a <- sum(tcrossprod(R, omega) * R) #quad.tdiag(omega, R) %>% sum()
    diagA <- colSums(crossprod(omega,fbar) * fbar) #quad.diag(omega, fbar) # diag( t(F) %*% omega %*% F)
    part_b <- t(t(varlbar)*diagA) %>% sum()
    ell <- -0.5*(part_a + part_b)
  }else{
    part_a <- map_dbl(seq(n), function(i){
      #quad.tform(omega[[i]], R[i,,drop = FALSE])
      r <- R[i,,drop = FALSE]
      tcrossprod(r, tcrossprod(r, omega[[i]]))
    }) %>% sum()

    diagA <- map(seq(n), function(i){
      colSums(crossprod(omega[[i]],fbar) * fbar)
    }) %>% unlist() %>% matrix(ncol = k, byrow = T)

    part_b <- sum(varlbar*diagA)
    #part_b <- t(t(varlbar)*diagA) %>% sum()
    ell <- -0.5*(part_a + part_b)
  }
  return(ell)
}

ell_const <- function(omega, n){
  s_equal <- check_equal_omega(omega)
  if(s_equal){
    p <- nrow(omega)
    c <- -0.5*(p*log(2*base::pi) + log(solve(det(omega))))
    return(n*as.numeric(c))
  }else{
    p <- nrow(omega[[1]])
    c <- sapply(omega, function(o){
      log(solve(det(o)))
    }) %>% sum()
    c <- -0.5*(c + n*p*log(2*base::pi))
    return(as.numeric(c))
  }
}

# computes
#  log p(yi | li, fbar, omega_i, mu_k, v_k) with
# p(yi | li fbar) ~ N(li t(F), omega^{-1})
# li ~ N(mu_k, diag(v_k))
ll_component <- function(Y, fbar, omega, mu_k, v_k){
  n <- nrow(Y)
  p <- ncol(Y)
  s_equal <- check_equal_omega(omega)
  check_matrix(fbar, "fbar", p, p)
  check_omega(omega, n, p, s_equal)
  stopifnot(length(mu_k) == p & length(v_k) == p)
  Vf <- fbar %*% diag(v_k) %*% t(fbar)
  fmu <- fbar %*% matrix(mu_k, nrow = p)
  R <- t(t(Y) - as.vector(fmu))
  if(s_equal){
    V <- Vf + solve(omega)
    Vinv <- solve(V)
    ll <- -0.5*rowSums(tcrossprod(R, Vinv) *R)#quad.tdiag(omega, R)
    c <- -0.5*(p*log(2*base::pi) + log(det(V)))
    ll <- ll + c
    #ll <- apply(Y, 1, function(y){ dmvnorm(y, mean  = fmu, sigma = V, log = TRUE)})
  }else{
    V <- lapply(omega, function(o){Vf + solve(o)})
    Vinv <- lapply(V, function(v){solve(v)})
    ll <-  map_dbl(seq(n), function(i){
      r <- R[i,,drop = FALSE]
      tcrossprod(r, tcrossprod(r, Vinv[[i]]))
    })
    c <-  -0.5*(p*log(2*base::pi) + sapply(V, function(v){log(det(v))}))
    ll <- ll + c
  }
  return(ll)
}


calc_ll <- function(Y, fbar, omega, ghat){
  n <- nrow(Y)
  p <- ncol(Y)
  ncomp <- lapply(ghat, function(g){which(g$pi > 0)})
  mean_zero <- sapply(ghat, function(g){all(g$mean == 0)}) %>% max
  if(!mean_zero) mu <- lapply(ghat, function(g){g$mean}) %>% expand.grid
  s <- lapply(ghat, function(g){g$sd[g$pi > 0]}) %>% expand.grid
  lpi <- lapply(ghat, function(g){log(g$pi[g$pi > 0])}) %>% expand.grid %>% rowSums
  nc <- nrow(s)
  if(mean_zero){
    ll <- sapply(seq(nc), function(i){
      ll_component(Y, fbar, omega, rep(0, p), s[i,]^2 )
    })
  }else{
    ll <- sapply(seq(nc), function(i){
      ll_component(Y, fbar, omega, mu[i,], s[i,]^2 )
    })
  }
  ll <- t(t(ll) + lpi)
  li <- apply(ll, 1, logSumExp)
  return(sum(li))
}

# computes
#  log p(yi | li, fbar, omega_i, mu_k, v_k) with
# p(yi | li fbar) ~ N(li t(F), omega^{-1})
# li ~ N(mu_k, L)
# L= diag(s_k)R diag(s_k)
ll_componentR <- function(Y, fbar, omega, mu_k, v_k, R){
  n <- nrow(Y)
  p <- ncol(Y)
  s_equal <- check_equal_omega(omega)
  check_matrix(fbar, "fbar", p, p)
  check_omega(omega, n, p, s_equal)
  stopifnot(length(mu_k) == p & length(v_k) == p)

  S <- diag(sqrt(v_k))
  L <- S %*% R %*% S

  Vf <- fbar %*% L %*% t(fbar)
  fmu <- fbar %*% matrix(mu_k, nrow = p)
  R <- t(t(Y) - as.vector(fmu))
  if(s_equal){
    V <- Vf + solve(omega)
    Vinv <- solve(V)
    ll <- -0.5*rowSums(tcrossprod(R, Vinv) *R)#quad.tdiag(omega, R)
    c <- -0.5*(p*log(2*base::pi) + log(det(V)))
    ll <- ll + c
    #ll <- apply(Y, 1, function(y){ dmvnorm(y, mean  = fmu, sigma = V, log = TRUE)})
  }else{
    V <- lapply(omega, function(o){Vf + solve(o)})
    Vinv <- lapply(V, function(v){solve(v)})
    ll <-  map_dbl(seq(n), function(i){
      r <- R[i,,drop = FALSE]
      tcrossprod(r, tcrossprod(r, Vinv[[i]]))
    })
    c <-  -0.5*(p*log(2*base::pi) + sapply(V, function(v){log(det(v))}))
    ll <- ll + c
  }
  return(ll)
}

calc_llR <- function(Y, fbar, omega, ghat, R){
  n <- nrow(Y)
  p <- ncol(Y)
  ncomp <- lapply(ghat, function(g){which(g$pi > 0)})
  mean_zero <- sapply(ghat, function(g){all(g$mean == 0)}) %>% max
  if(!mean_zero) mu <- lapply(ghat, function(g){g$mean}) %>% expand.grid
  s <- lapply(ghat, function(g){g$sd[g$pi > 0]}) %>% expand.grid
  lpi <- lapply(ghat, function(g){log(g$pi[g$pi > 0])}) %>% expand.grid %>% rowSums
  nc <- nrow(s)
  if(mean_zero){
    ll <- sapply(seq(nc), function(i){
      ll_componentR(Y, fbar, omega, rep(0, p), s[i,]^2, R )
    })
  }else{
    ll <- sapply(seq(nc), function(i){
      ll_componentR(Y, fbar, omega, mu[i,], s[i,]^2, R )
    })
  }
  ll <- t(t(ll) + lpi)
  li <- apply(ll, 1, logSumExp)
  return(sum(li))
}
