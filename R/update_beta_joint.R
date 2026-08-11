
#'@export
update_beta_joint <- function(dat,
                              j=1,
                              ix = NULL,
                              ii = NULL,
                              return_W = FALSE,
                              cond_num = 1e10){

  # j is the row of F that we will update. F is p rows by k columns
  # for factors version, abar and lbar are both n by k, lbar and abar are the same

  p <- dat$p
  k <- dat$k
  n <- dat$n

  prior_precision <- dat$beta$prior_precision

  if(is.null(ix)){
    ix <- seq(k)[-j]
  }else{
    stopifnot(all(ix %in% seq(k)))
    stopifnot(!any(duplicated(ix)))
    #ix <- sort(ix)
  }
  m <- length(ix)
  if(is.null(prior_precision)){
    T0 <- matrix(0, nrow = m, ncol = m)
  }else{
    T0 <- prior_precision[ii,ii, drop = FALSE]
  }
  Va <- dat$l$a2bar - (dat$l$abar^2)

  if(dat$s_equal){
    A <- t(dat$l$abar) %*% dat$l$abar + diag(colSums(Va))
    Astar <- dat$G %*% A %*% t(dat$G)

    Rfull <- dat$omega[j,j]*Astar  # W in the manuscript k by k
    a10 <- colSums(dat$l$lbar *rowSums(t(t(dat$Y)*dat$omega[,j]))) # length k
    a20 <- lapply(seq(p)[-j], function(jj){
      # (k by k ) %*% (k by 1)*(p by 1)
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
    R <- Rfull[ix,ix, drop = FALSE]
    R12 <- Rfull[ix,-ix, drop = FALSE]
    a <- afull[ix] - R12 %*% t(dat$f$fbar[j,-ix,drop = FALSE])
  }else{
    R <- Rfull
    a <- afull
  }

  remove_suggest <- NULL
  evR <- eigen(R, only.values = TRUE)$values
  condR <- abs(max(evR)/min(evR))
  if(any(evR < 0) | condR > cond_num){
    info_abar <- colSums(dat$l$a2bar)
    if(dat$is_nesmr){
      worst_abar <- which.min(info_abar)
      remove_suggest <- which.max(abs(dat$G[,worst_abar]))
      warning("There is not enough independent genetic information to estimate all trait effects.
              This will result in some very large standard errors.  I recommend removing trait ", remove_suggest, ".\n")
    }else if(dat$is_factors){ ## Factors case
      worst_abar <- which.min(info_abar[-1]) + 1
      remove_suggest <- which.max(abs(dat$G[,worst_abar]))
      if(remove_suggest == 2){
        warning("There may not be enough independent genetic information about your primary exposure. This could happen if there are
                very few or very weak instruments for X. This will result in very large standard errors.")
      }else{
        warning("There is not enough independent genetic information to estimate all factor effects.
              I am going to remove factor ", remove_suggest-2, ".\n")
      }
    }else{ # MVMR case
      info_lbar <- colSums(dat$l$l2bar)
      remove_suggest <- which.min(info_lbar[-1]) + 1 ## do not check Y
      warning("There is not enough independent genetic information to estimate all trait effects.
              This will result in some very large standard errors.  I recommend removing exposure trait ", remove_suggest-1, ".\n")

    }
    warning("Projecting internal R to nearest PD matrix in beta update.\n")
    R <- Matrix::nearPD(R, posd.tol = 1/cond_num)$mat
  }

  S <- solve(R + T0)
  mu <- S %*% a

  if(return_W){
    return(list(m = mu, S = S, W = R, b = a, remove_suggest = remove_suggest))
  }
  return(list(m = mu, S = S, remove_suggest = remove_suggest))
}

## This only gets used for nesmr and factor problems
update_beta_full_joint <- function(dat, cond_num = 1e10){

  p <- dat$p
  n <- dat$n
  k <- dat$k

  ix <- p*(dat$beta$beta_k-1) + dat$beta$beta_j
  ix <- ix[!dat$beta$fix_beta]
  m <- length(ix)

  prior_precision <- dat$beta$prior_precision
  if(is.null(prior_precision)){
    T0 <- matrix(0, nrow = m, ncol = m)
  }else{
    T0 <- prior_precision
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
    # Rfull <- lapply(seq(n), function(i){
    #   l <- dat$l$abar[i,]
    #   a <- outer(l, l) + diag(Va[i,], nrow = k)
    #   kronecker(tcrossprod(dat$G, tcrossprod(dat$G, a)), dat$omega[[i]]) ## kronecker(G %*% a %*% t(G), O)
    # }) %>% Reduce(`+`, .)
    Gt <- t(dat$G)
    nG <- nrow(dat$G)
    nO <- nrow(dat$omega[[1]])
    Rfull <- matrix(0, nG * nO, nG * nO)
    for (i in seq(n)) {
      l <- dat$l$abar[i,]
      a <- tcrossprod(l) + diag(Va[i,], nrow = k)
      A <- dat$G %*% tcrossprod(a, Gt)
      O <- dat$omega[[i]]
      for (r in seq_len(nG)) {
        ri <- ((r-1)*nO + 1):(r*nO)
        for (s in seq_len(nG)) {
          si <- ((s-1)*nO + 1):(s*nO)
          Rfull[ri, si] <- Rfull[ri, si] + A[r, s] * O
        }
      }
    }
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

  remove_suggest <- NULL
  evR <- eigen(R, only.values = TRUE)$values
  condR <- abs(max(evR)/min(evR))
  if(any(evR < 0) | condR > cond_num){
    info_abar <- colSums(dat$l$abar^2)
    if(dat$is_nesmr){
      worst_abar <- which.min(info_abar)
      remove_suggest <- which.max(abs(dat$G[,worst_abar]))
      warning("There is not enough independent genetic information to estimate all trait effects.
              This will result in some very large standard errors.  I recommend removing trait ", remove_suggest, ".\n")
    }else if(dat$is_factors){ ## Factors case
      worst_abar <- which.min(info_abar[-1]) + 1
      remove_suggest <- which.max(abs(dat$G[,worst_abar]))
      if(remove_suggest == 2){
        warning("There may not be enough independent genetic information about your primary exposure. This could happen if there are
                very few or very weak instruments for X. This will result in very large standard errors.")
      }else{
        warning("There is not enough independent genetic information to estimate all factor effects.
              I am going to remove factor ", remove_suggest-2, ".\n")
      }
    }else{ # MVMR case
      info_lbar <- colSums(dat$l$l2bar)
      remove_suggest <- which.min(info_lbar[-1]) + 1 ## do not check Y
      warning("There is not enough independent genetic information to estimate all trait effects.
              This will result in some very large standard errors.  I recommend removing exposure trait ", remove_suggest-1, ".\n")

    }

    warning("Projecting internal R to nearest PD matrix in beta update.\n")
    R <- Matrix::nearPD(R, posd.tol = 1/cond_num)$mat
  }
  S <- solve(R + T0)
  mu <- S %*% a
  return(list(m = mu, S = S, remove_suggest = remove_suggest))
}
