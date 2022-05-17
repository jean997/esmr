
#'@export
eb_mr <- function(beta_hat_Y, se_Y, beta_hat_X, se_X,
                  R = NULL,
                  ebnm_fn = flashier:::as.ebnm.fn(prior_family = "point_normal", optmethod = "nlm"),
                  max_iter = 100,
                  seed = 123, sigma_beta = Inf,
                  tol = 1e-7,
                  lfsr_thresh = 1, pval_thresh =1, pval_select,
                  beta_m_init, which_beta,
                  fix_beta = FALSE,
                  beta_joint = FALSE){

  set.seed(seed)

  n <- length(beta_hat_Y)
  stopifnot(length(se_Y) == n)

  if("matrix" %in% class(beta_hat_X)){
    stopifnot(nrow(beta_hat_X) == n)
    stopifnot(nrow(se_X) == n)
    p <- ncol(beta_hat_X) + 1
  }else{
    stopifnot(length(beta_hat_X) == n)
    stopifnot(length(se_X) == n)
    p <- 2
  }
  if(missing(R) | is.null(R)){
    R <- diag(p)
  }

  if(pval_thresh < 1){
    if(missing(pval_select)) stop("If pval_thresh is <1 please provide pval_select")
    ix_pval <- which(pval_select < pval_thresh)
  }else{
    ix_pval <- 1:n
  }

  if(missing(which_beta) | beta_joint){
    beta_j <- rep(1, p-1)
    beta_k <- 2:p
    if(missing(beta_m_init)){
      beta_m <- rep(0, p-1)
      beta_s <- rep(0, p-1)
    }else{
      stopifnot(length(beta_m_init) == p-1)
      beta_m <- beta_m_init
      beta_s <- rep(0, p-1)
    }
  }else if(!beta_joint){
    beta_j <- which_beta[,1]
    beta_k <- which_beta[,2]
    if(missing(beta_m_init)){
      beta_m <- rep(0, length(beta_j))
      beta_s <- rep(0, length(beta_j))
    }else{
      stopifnot(length(beta_m_init) == length(beta_j))
      beta_m <- beta_m_init
      beta_s <- rep(0, length(beta_j))
    }
  }



  Y <- cbind(beta_hat_Y, beta_hat_X)
  S <- cbind(se_Y, se_X)

  s_equal <- apply(S, 2, function(x){all(x == x[1])}) %>% all()
  if(s_equal){
    s <- S[1,]
    omega <- solve(diag(s) %*% R %*% diag(s))
  }else{
    omega <- apply(S, 1, function(s){
      solve(diag(s) %*% R %*% diag(s))
    }, simplify = FALSE)
  }

  lbar <- l2bar <- matrix(0, nrow = n, ncol = p)
  # l2bar <- l_init %>% map(function(l){l$posterior$second_moment}) %>% do.call(cbind, .)
  i <- 1

  make_fbar <- function(beta_m, beta_s, beta_j, beta_k, p){
    fbar <- f2bar <- diag(p)
    for(i in seq(length(beta_j))) fbar[beta_j[i],beta_k[i]] <- beta_m[i]
    for(i in seq(length(beta_j))) f2bar[beta_j[i],beta_k[i]] <- beta_m[i]^2 + beta_s[i]^2
    return(list(fbar = fbar, f2bar = f2bar))
  }

  ff <- make_fbar(beta_m, beta_s, beta_j, beta_k, p)
  fbar <- ff$fbar
  f2bar <-  ff$f2bar

  check <- 1
  kl <- c()
  ll <- c()
  obj <- c()
  obj_old <- -Inf
  l_update <- list()
  while(i < max_iter & check > tol){
    # if(i == 1) kk <- 1
    #   else kk <- 1
    # alpha and gamma updates
    #for(kk in 1:kk){

    for(j in seq(p)){
        R_j <- Y - (lbar[,-j,drop=FALSE] %*% t(fbar[,-j,drop=FALSE]))
        lu <- update_l_k(R_j, fbar[,j], f2bar[,j], omega, ebnm_fn)
        lbar[,j] <- lu$posterior$mean
        l2bar[,j] <- lu$posterior$second_moment
        l_update[[j]] <- lu
        #kl <- c(kl, map(l_update, "KL") %>% unlist() %>% sum())
        #ll <- c(ll, calc_ll(Y, lbar, l2bar, beta_m, beta_s, omega))
    }
      #lbar <- l_update %>% map(function(l){l$posterior$mean}) %>% do.call(cbind, .)
      #l2bar <- l_update %>% map(function(l){l$posterior$second_moment}) %>% do.call(cbind, .)
    lfsr <- l_update %>% map(function(l){l$posterior$lfsr}) %>% do.call(cbind, .)
    l_ghat <- l_update %>% map(function(l){l$fitted_g})



    # beta update
    if(!fix_beta & !beta_joint){
      R_k <- map(seq(p), function(k){
        Y[ix_pval,] - (lbar[ix_pval,-k,drop=FALSE]%*%t(fbar[,-k,drop=FALSE]))
      })
      if(s_equal){
        beta_upd <- map(seq(length(beta_j)), function(r){
          k <- beta_k[r]
          j <- beta_j[r]
          ixk <- which(lfsr[ix_pval,k] <= lfsr_thresh)
          update_beta_k(R_k = R_k[[k]][ixk,], j=j, k=k,
                        lbar=lbar[ix_pval[ixk],], l2bar=l2bar[ix_pval[ixk],],
                        omega = omega, fbar = fbar,
                        sigma_beta = sigma_beta)
        })
      }else{
        beta_upd <- map(seq(length(beta_j)), function(r){
          k <- beta_k[r]
          j <- beta_j[r]
          ixk <- which(lfsr[,k] < lfsr_thresh)
          update_beta_k(R_k = R_k[[k]][ixk,], j=j, k=k, lbar=lbar[ixk,], l2bar=l2bar[ixk,],
                        omega = omega[ixk], fbar = fbar,
                        sigma_beta = sigma_beta)
        })
      }
      beta_m <- map(beta_upd, "m") %>% unlist()
      beta_s <- map(beta_upd, "s") %>% unlist()
      ff <- make_fbar(beta_m, beta_s, beta_j, beta_k, p)
      fbar <- ff$fbar
      f2bar <-  ff$f2bar
    }else if(!fix_beta & beta_joint){
      beta_upd <- update_beta_joint(Y[ix_pval,], lbar[ix_pval,], l2bar[ix_pval,], omega)
      beta_m <- beta_upd$m
      beta_s <- sqrt(diag(beta_upd$S))
      beta_var <- beta_upd$S
      ff <- make_fbar(beta_m, beta_s, beta_j, beta_k, p)
      fbar <- ff$fbar
      f2bar <-  ff$f2bar
    }

    kl <- c(kl, map(l_update, "KL") %>% unlist() %>% sum())
    #ll <- c(ll, calc_ll(Y, lbar, l2bar, beta_m, beta_s, omega))
    #obj <- ll + kl
    obj <- kl
    obj_new <- obj[length(obj)]
    check <- obj_new - obj_old
    obj_old <- obj_new

    #if(check < -1e-5) warning("Objective decreased on iteration ", i, ".\n")
    check <- abs(check)
    cat(i, ": ", obj_new, " ", beta_m, " ", beta_s, "\n")
    i <- i + 1
  }

  result <- list(beta_m = beta_m, beta_s = beta_s, Y = Y, S = S,
                 lbar = lbar, l2bar = l2bar, fbar = fbar,
                 l_fit = l_update, omega = omega, obj = obj, ll = ll, kl = kl)
  if(beta_joint) result$beta_var <- beta_var
  return(result)
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
