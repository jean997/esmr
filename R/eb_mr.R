
#'@export
eb_mr <- function(beta_hat_Y, se_Y, beta_hat_X, se_X,
                  R = NULL,
                  ebnm_fn = flashier:::as.ebnm.fn(prior_family = "point_normal", optmethod = "nlm"),
                  max_iter = 100,
                  seed = 123, sigma_beta = Inf,
                  tol = 1e-7,
                  lfsr_thresh = 1, pval_thresh =1, pval_select,
                  beta_m_init, which_beta,
                  fix_beta = FALSE){

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

  if(missing(which_beta)){
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
  }else{
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
  obj <- -Inf

  while(i < max_iter & check > tol){
    if(i == 1) kk <- 3
      else kk <- 1
    # alpha and gamma updates
    for(kk in 1:kk){
      l_update <- map(rev(seq(p)), function(j){
        R_j <- Y - (lbar[,-j,drop=FALSE] %*% t(fbar[,-j,drop=FALSE]))
        update_l_k(R_j, fbar[,j], f2bar[,j], omega, ebnm_fn)
      }) %>% rev()
      lbar <- l_update %>% map(function(l){l$posterior$mean}) %>% do.call(cbind, .)
      l2bar <- l_update %>% map(function(l){l$posterior$second_moment}) %>% do.call(cbind, .)
      lfsr <- l_update %>% map(function(l){l$posterior$lfsr}) %>% do.call(cbind, .)
      l_ghat <- l_update %>% map(function(l){l$fitted_g})
    }

    R_k <- map(seq(p), function(k){
      Y[ix_pval,] - (lbar[ix_pval,-k,drop=FALSE]%*%t(fbar[,-k,drop=FALSE]))
    })
    # beta update
    if(!fix_beta){
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
    }

    if(!is.null(ebnm_fn)){
      obj_old <- obj
      obj <- map(l_update, "KL") %>% unlist() %>% sum()
      check <- obj - obj_old
      #if(check < 0) warning("Objective increased on iteration ", i, ".\n")
      check <- abs(check)
    }else{
      obj <- NA
      check <- 1
    }
    cat(i, ": ", obj, " ", beta_m, " ", beta_s, "\n")
    i <- i + 1
  }
  fitted_values <- lbar %*% t(fbar)
  result <- list(beta_m = beta_m, beta_s = beta_s, Y = Y, S = S,
                 lbar = lbar, fbar = fbar,
                 l_fit = l_update, fitted_values = fitted_values)
  return(result)
}


calc_ll <- function(Y, lbar, l2bar, fbar, f2bar, omega, l_ghat){
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
