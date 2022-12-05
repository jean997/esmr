#'@export
update_beta_k_future <- function(j, k, dat){
  s_equal <- check_equal_omega(dat$omega)

  Vl <- dat$l$l2bar_o - (dat$l$lbar_o^2)

  if(s_equal){
    A <- t(dat$l$lbar_o) %*% dat$l$lbar_o + diag(colSums(Vl))
    Astar <- dat$G %*% A %*% t(dat$G)

    s2inv <- (1/dat$sigma_beta^2) + Astar[k,k]*dat$omega[j,j]
    s2 <- 1/s2inv

    m_pt1 <- sum(dat$l$lbar[,k]*t(dat$omega[j,,drop =FALSE]%*%t(dat$Y)))
    m_pt2 <- sum(Astar[k, -k] * (dat$omega[j,] %*% dat$f$fbar_o[,-k]))
    m_pt3 <- Astar[k,k]*sum(dat$f$fbar_o[-j,k]*dat$omega[j,-j])

    m <- (m_pt1 - m_pt2 - m_pt3) * s2
  }else{
    stop("Not implemented yet.")
    ojj <- map(dat$omega, function(o){o[j,j]}) %>% unlist()

    s2inv <- (1/dat$sigma_beta^2) + sum(dat$l$l2bar[,k]*ojj)
    s2 <- 1/s2inv

    Oj <- map(omega, function(o){o[j,]}) %>% unlist() %>% matrix(nrow = n, byrow = TRUE)
    Ojmj <- map(omega, function(o){o[j,-j]}) %>% unlist() %>% matrix(nrow = n, byrow = TRUE)
    m_over_s2_pt1 <- sum(lbar[,k]*diag(Oj%*%t(R_k)))

    m_over_s2_pt2 <- sum(l2bar[,k]*t(t(fbar[-j,k,drop=FALSE])%*%t(Ojmj)))

    m <- (m_over_s2_pt1 - m_over_s2_pt2) * s2
  }

  return(list(s = sqrt(s2), m = m))

}

update_beta_sequential_future <- function(dat, jj, fix_beta = NULL){
  p <- dat$p

  beta_j <- dat$beta$beta_j
  beta_k <- dat$beta$beta_k

  fbar_o <- dat$f$fbar_o
  f2bar_o <- dat$f$f2bar_o
  beta_m <- dat$beta$beta_m
  beta_s <- dat$beta$beta_s

  if(missing(jj)){
    coords <- seq(length(beta_j))
  }else{
    coords <- jj
  }
  if(!is.null(fix_beta)){
    if(length(fix_beta) == 1){
      fix_beta <- rep(fix_beta, length(coords))
    }else if(length(fix_beta) != length(coords)){
      stop("fix_beta and expected to have length 1 or ", length(coords), ". Found ", length(fix_beta), ".\n")
    }
    stopifnot(class(fix_beta) == "logical")
    coords <- coords[!fix_beta]
  }

  for(i in coords){
    k <- beta_k[i]
    j <- beta_j[i]

    #R_k <- dat$Y - dat$l$lbar[,-k,drop=FALSE]%*%t(fbar_o[,-k,drop=FALSE])
    #b <- update_beta_k(R_k = R_k, j=j, k=k,
    #                   lbar=dat$l$lbar, l2bar=dat$l$l2bar,
    #                   omega = dat$omega, fbar = fbar_o,
    #
    b <- update_beta_k_future(j,k,dat)
    beta_m[i] <- b$m
    beta_s[i] <- b$s
    f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
    fbar_o <- f$fbar_o
    f2bar_o <- f$f2bar_o
  }
  dat$beta$beta_m <- beta_m
  dat$beta$beta_s <- beta_s
  dat$f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
  return(dat)
}

#'@export
update_beta_joint_future <- function(dat, j=1, ix = NULL, prior_cov = NULL){
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
    stopifnot(any(duplicated(ix)))
    ix <- sort(ix)
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
