
one_elbo <- function(beta, pi0, s1, abari, a2bari, Yi, omegai, G, fbar, ebnm_fn, j = 1, ix = NULL){
  p <- ncol(Yi)
  k <- ncol(G)
  if(is.null(ix)){
    ix <- seq(p)[-j]
  }
  fbar[j, ix] <- beta
  fgbar <- fbar %*% G

  ## step 3 calc_ll2
  ll <- calc_ell2(Yi, abari, a2bari, fgbar, omegai)
  #return(ll)

  ghat <- lapply(1:k, function(kk){
    ashr::normalmix(pi = c(pi0[kk], 1-pi0[kk]), mean = c(0, 0), sd = c(0, s1[kk]))
  })
  ## calculate KL
  kl <- sapply(1:k, function(kk){
   a <- 1/ghat[[kk]]$sd[2]^2
   w <- ghat[[kk]]$pi[2]
   x_s <- get_x_k(Yi, abari, fgbar, kk , omegai, ebnm_fn)
   calc_kl_k(x_s$x, x_s$s, w, a, 0, abari[,kk], a2bari[,kk] )
  }) |> sum()
  return(ll + kl)
}

one_elbo_gradA <- function(beta, dat, ix){
  w <- update_beta_joint(dat = subset_data(dat, ix), return_W = T)
  return(w$b - w$W %*% matrix(beta, ncol = 1))
}

one_elbo_hessA <- function(beta, dat, ix){
  w <- update_beta_joint(dat = subset_data(dat, ix), return_W = T)
  return(-w$W)
}

one_elboA <- function(beta, pi0, s1, abari, a2bari, Yi, omegai, G, fbar, ebnm_fn, j = 1, ix = NULL){
  p <- col(Yi)
  k <- ncol(G)
  if(is.null(ix)){
    ix <- seq(p)[-j]
  }
  fbar[j, ix] <- beta
  fgbar <- fbar %*% G


  ## step 3 calc_ll2
  ll <- calc_ell2(Yi, abari, a2bari, fgbar, omegai)
  return(ll)
}

one_elboB <- function(beta, pi0, s1, abari, a2bari, Yi, omegai, G, fbar, ebnm_fn, j = 1, ix = NULL){
  p <- ncol(Yi)
  k <- ncol(G)
  if(is.null(ix)){
    ix <- seq(p)[-j]
  }
  fbar[j, ix] <- beta
  fgbar <- fbar %*% G

  ghat <- lapply(1:k, function(kk){
    ashr::normalmix(pi = c(pi0[kk], 1-pi0[kk]), mean = c(0, 0), sd = c(0, s1[kk]))
  })
  ## calculate KL
  kl <- sapply(1:k, function(kk){
    a <- 1/ghat[[kk]]$sd[2]^2
    w <- ghat[[kk]]$pi[2]
    x_s <- get_x_k(Yi, abari, fgbar, kk , omegai, ebnm_fn)
    calc_kl_k(x_s$x, x_s$s, w, a, 0, abari[,kk], a2bari[,kk] )
  }) |> sum()
  return(kl)
}


get_x_k <- function(Y, abar, fgbar, k , omega, ebnm_fn){
  R_k <- Y - (abar[,-k,drop=FALSE] %*% t(fgbar[,-k,drop=FALSE]))
  x_s <- update_l_k(R_k, fgbar[,k], fgbar[,k]^2, omega, ebnm_fn,
                    return_x_s = TRUE,
                    g_init = ghat[[j]], fix_g = TRUE)
  return(x_s)
}


sandwich_cov <- function(beta, pi0, s1, j = 1, dat,
                         subset_ix = NULL,
                         verbose=TRUE, nodes=1) {

  if(!is.null(subset_ix)){
    dat <- subset_data(dat, subset_ix)
  }

  if(missing(beta)){
    beta <- dat$beta$beta_m
  }
  if(missing(pi0)){
    pi0 <- sapply(dat$l$g_hat, function(x){x$pi[1]})
  }
  if(missing(s1)){
    s1 <- sapply(dat$l$g_hat, function(x){x$sd[2]})
  }

  theta <- c(beta, log(pi0/(1-pi0)), -log(1/s1^2))

  k <- ncol(dat$G)
  psi_ln <- 2*k
  N <- nrow(dat$Y)
  p <- ncol(dat$Y)
  ix1 <- seq(length(beta))
  ix1_ln <- length(beta)
  ix2 <- seq(4*k) + length(beta)






  #abs <- lapply(1:N, function(i) {
  abs <- mclapply(1:N, function(i) {
    if(verbose ) cat(i, " of ", N, " \n")


    grad <- one_elbo_gradA(beta, dat, i)

    hessA <- one_elbo_hessA(beta, dat, i)

    va <- dat$l$a2bar[i,] - (dat$l$abar[i,]^2)
    g_a <- -log(1/va)
    thetapsi <-  c(theta, dat$l$abar[i,], g_a)


    fB <- function(thetapsi){
      beta <- thetapsi[ix1]
      alpha <- thetapsi[p:(p + k -1)]
      pi0  <- 1 / (exp(-alpha) + 1) # pi0 <- expit(alpha), alpha <- logit(pi0)
      gamma <- thetapsi[(p + k):(p + 2*k-1)]
      # gamma <- -log( 1/s1^2), s1 <- sqrt(1/exp(-gamma))
      s1 <- sqrt(1/exp(-gamma))
      abari <- matrix(thetapsi[(p + 2*k):(p + 3*k -1 )], nrow =1 )
      gamma_a <- thetapsi[(p + 3*k):(p + 4*k -1)]
      va <- 1/exp(-gamma_a)
      a2bari <- va + abari^2
      a2bari <- matrix(a2bari, nrow = 1)
      one_elboB(beta, pi0, s1, abari, a2bari, Yi = dat$Y[i,,drop=F],
                omegai = dat$omega[[i]], G = dat$G, fbar = dat$f$fbar, ebnm_fn = dat$ebnm_fn)
    }


    #g2 <- maxLik::numericGradient(f, t0 = thetapsi )
    hessB <- maxLik::numericHessian(fB,t0 =thetapsi)
    hess <- hessB
    hess[ix1, ix1] <- hessB[ix1, ix1] + hessA

    X <- hess[ix2, ix2]
    eX <- eigen(X)
    sX <- with(eX, vectors %*% diag(1/values) %*% t(vectors))
    #sX <- solve(X)
    Ai <- hess[ix1, ix1] - hess[ix1, ix2] %*% sX %*% hess[ix2, ix1]

    Bi <- outer(grad, grad)
    c(c(Ai), c(Bi))
  }, mc.cores=nodes)
  abs <- matrix(unlist(abs), nrow=length(abs), byrow=TRUE)
  abhat <- colMeans(abs)
  ahat <- matrix(abhat[1:ix1_ln^2], nrow=ix1_ln)
  bhat <- matrix(abhat[-(1:ix1_ln^2)], nrow=ix1_ln)
  return(list(ahat=ahat, bhat=bhat, sand_cov=solve(ahat) %*% bhat %*% solve(ahat) / N))
}


get_abar_a2bar <- function(Y, omega, G, ghat, fbar, ebnm_fn, tol = 1e-8, nmax = 100){
  n <- nrow(Y)
  k <- ncol(G)
  abar <- matrix(0, nrow = n, ncol = k)
  a2bar <- matrix(0, nrow = n, ncol = k)
  fgbar <- fbar %*% G
  kl <- rep(0, k)

  abar0 <- a2bar0 <- abar
  done <- FALSE
  i <- 0
  check <- Inf
  while(!done){
    i <- i + 1
    #cat(i, " ", check, "\n")
    for(j in 1:k){
      R_j <- Y - (abar[,-j,drop=FALSE] %*% t(fgbar[,-j,drop=FALSE]))
      lu <- update_l_k(R_j, fgbar[,j], fgbar[,j]^2, omega, ebnm_fn, g_init = ghat[[j]], fix_g = TRUE)
      abar[lu$posterior$index,j] <- lu$posterior$mean
      a2bar[lu$posterior$index,j] <- lu$posterior$second_moment
      kl[j] <- lu$KL
    }
    check1 <- max(abs(abar0 - abar))
    check2 <- max(abs(a2bar0 - a2bar))
    check <- max(check1, check2)
    if(check < tol) done <- TRUE
    if(i >= nmax ) done <- TRUE
    abar0 <- abar
    a2bar0 <- a2bar
  }
  return(list(abar = abar, a2bar = a2bar, kl = sum(kl)))
}
