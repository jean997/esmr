
#'@export
eb_mr <- function(beta_hat_Y, se_Y, beta_hat_X, se_X,
                  R = NULL,
                  ebnm_fn = flashier:::as.ebnm.fn(prior_family = "point_normal", optmethod = "nlm"),
                  max_iter = 100,
                  seed = 123, sigma_beta = Inf,
                  tol = default_precision(c(ncol(beta_hat_X)+1, nrow(beta_hat_X))),
                  pval_thresh =1, lfsr_thresh = 1,
                  beta_m_init = NULL, which_beta = NULL,
                  fix_beta = FALSE,
                  beta_joint = TRUE){

  set.seed(seed)

  dat <- set_data(beta_hat_Y, se_Y, beta_hat_X, se_X, R)
  n <- dat$n
  p <- dat$p

  dat$beta <- init_beta(p, which_beta, beta_joint, beta_m_init)
  dat$f_fun <- make_f_fun(p, dat$beta$beta_j, dat$beta$beta_k)

  dat$l <- init_l(n, p)
  dat$f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
  dat$fix_beta <- fix_beta
  dat$beta_joint <- beta_joint
  dat$ebnm_fn <- ebnm_fn
  dat$sigma_beta <- sigma_beta

  dat$pval <- with(dat, 2*pnorm(-abs(Y/S)))
  dat$pval <- apply(dat$pval[,-1], 1, min)
  dat$ix <- which(dat$pval < pval_thresh)

  dat$lfsr_thresh <- lfsr_thresh
  dat$pval_thresh <- pval_thresh

  dat <- ebmr_solve(dat, max_iter, tol )

  return(dat)
}


calc_ll <- function(Y, lbar, l2bar, beta_m, beta_s, omega){
  n <- nrow(Y)
  p <- ncol(Y)
  stopifnot(nrow(lbar) == n & ncol(lbar) == p)
  stopifnot(nrow(l2bar) == n & ncol(l2bar) == p)
  stopifnot(length(beta_m) == p-1 & length(beta_s) == p-1)


  #First col of lbar is alpha the rest are gammas
  if("matrix" %in% class(omega)){
    s_equal <- TRUE
  }else{
    stopifnot(class(omega) != "list")
    stopifnot(length(omega) == n)
    s_equal <- FALSE
  }
  beta_m <- as.vector(beta_m)
  b2 <- beta_m^2 + (beta_s)^2
  s2l <- l2bar - (lbar^2)
  r <- Y - lbar
  if(s_equal){
    beta_lbar <- t(t(lbar[,-1, drop = FALSE])*beta_m)
    G <- t(t(r)*omega[1,])
    gi <- rowSums(G)
    t1 <- sum(beta_lbar^2)*omega[1,1]
    t2 <- sum(t(s2l[,-1])*b2)*omega[1,1]
    t3 <- -2*sum(rowSums(beta_lbar)*gi)
    t4 <- 2*sum(t(s2l[,-1])*(beta_m*omega[1,-1]))
    t5 <- sapply(seq(n), function(i){
      r[i,,drop = FALSE]%*%omega%*%t(r[i,,drop=FALSE])
    }) %>% unlist() %>% sum()
    t6 <- sum(t(s2l)*diag(omega))
    t7 <- 0.5*sum(log(s2l))
    x <- t1 + t2 + t3 + t4 + t5 + t6
    return(-0.5*x)
  }
}

calc_ll0 <- function(Y, lbar, l2bar, fbar, f2bar, omega, l_ghat){
  n <- nrow(Y)
  p <- ncol(Y)
  stopifnot(nrow(lbar) == n & ncol(lbar) == p)
  stopifnot(nrow(l2bar) == n & ncol(l2bar) == p)
  stopifnot(ncol(fbar) ==p)


  #First col of lbar is alpha the rest are gammas
  if("matrix" %in% class(omega)){
    s_equal <- TRUE
  }else{
    stopifnot(class(omega) != "list")
    stopifnot(length(omega) == n)
    s_equal <- FALSE
  }

  Ybar <- lbar %*% t(fbar)

  Y2bar <- l2bar %*% t(f2bar)
  R <- Y - Ybar


  # Calculate E(log p(Y | l, f, omega))
  Fbar <- map(seq(p), function(k){
    Fb <- fbar[,k] %*% t(fbar[,k])
    diag(Fb) <- f2bar[,k]
    return(Fb)
  })

  L2 <- lbar %*% t(lbar)
  diag(L2) <- 0
  if(s_equal){
    ll_t1 <- map(seq(n), function(i){
      Y[i,,drop = FALSE] %*% omega %*% t(R[i,,drop = FALSE])
    }) %>% unlist() %>% sum()


    ll_t2 <- map(seq(n), function(i){
      Y2i <- Ybar[i,,drop = FALSE]%*%omega %*% t(Ybar[i,,drop = FALSE])
      sum(Y2i)
    }) %>% unlist() %>% sum()

    ll_t0 <- map(seq(n), function(i){
      Y2i <- Y[i,,drop = FALSE]%*%omega %*% t(Y[i,,drop = FALSE])
      sum(Y2i)
    }) %>% unlist() %>% sum()

  }else{

  }

  #log likelihood minus constant depending on omega
  ll <- -0.5*(ll_t0 -2*ll_t1 + ll_t2)

  gll <- map(seq(p), function(k){
    llk <- ashr:::lik_normalmix(pilik = l_ghat[[k]]$pi,
                                sdlik = l_ghat[[k]]$sd)
    sum(llk$lpdfFUN(lbar[,k]))
  }) %>% unlist() %>% sum()

  return(ll + gll)

}
