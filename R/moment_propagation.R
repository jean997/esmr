moment_propagation <- function(dat, nsamps = 1000){

  dat <- update_l_sequential(dat, seq(dat$p), dat$g_init, dat$fix_g, return_sampler = TRUE)

  asamps <- lapply(seq(dat$p), function(i){
    dat$l$sampler[[i]](nsamps)
  })
  nb <- length(dat$beta$beta_j)
  dat$beta$V <- matrix(0, nrow = nb, ncol = nb)

  if(dat$R_is_id | length(unique(dat$beta$beta_j)) == 1){
    # if all omega are diagonal or only estimating one row, update F by rows
    jj <- unique(dat$beta$beta_j)
    for(j in jj){
      ii <- which(dat$beta$beta_j == j & !dat$beta$fix_beta)
      if(length(ii) == 0) next
      ix <- dat$beta$beta_k[ii]
      beta_upd <- update_beta_joint_mp(dat, asamps, j = j, ix = ix, ii = ii)

      dat$beta$beta_m[ii] <- beta_upd$m
      dat$beta$beta_s[ii] <- sqrt(diag(beta_upd$S))
      dat$beta$V[ii,ii] <- beta_upd$S

      dat$f <- make_f(dat)
    }
  }else{
    e_ix <- which(!dat$beta$fix_beta)
    ub <- update_beta_full_joint_mp(dat, asamps)
    dat$beta$beta_m[e_ix] <- ub$m
    dat$beta$V[e_ix,e_ix] <- ub$S
    dat$beta$beta_s[e_ix] <- sqrt(diag(ub$S))

    dat$f <- make_f(dat)
  }
}


#'@export
update_beta_joint_mp <- function(dat, asamps,
                                 j=1, ix = NULL, ii = NULL,
                                return_W = FALSE){
  p <- dat$p
  n <- dat$n
  ns <- nrow(asamps[[1]])
  prior_precision <- dat$beta$prior_precision

  if(is.null(ix)){
    ix <- seq(p)[-j]
  }else{
    stopifnot(all(ix %in% seq(p)))
    stopifnot(!any(duplicated(ix)))
    #ix <- sort(ix)
  }
  m <- length(ix)
  if(is.null(prior_precision)){
    T0 <- matrix(0, nrow = m, ncol = m)
  }else{
    # T0 <- check_matrix(prior_precision, m, m)
    T0 <- prior_precision[ii,ii]
  }


  if(dat$s_equal){
    samps <- lapply(seq(ns), function(i){
      asamp <- map(asamps, function(a){a[i,]}) %>% do.call(cbind, .)
      lsamp <- asamp %*% t(dat$G)
      A <- t(asamp) %*% asamp
      Astar <- dat$G %*% A %*% t(dat$G)

      Rfull <- dat$omega[j,j]*Astar  # W in the manuscript
      a10 <- colSums(lsamp *rowSums(t(t(dat$Y)*dat$omega[,j])))
      a20 <- lapply(seq(p)[-j], function(jj){
        Astar%*% t(dat$f$fbar[jj,,drop = FALSE])*dat$omega[j,jj]
      }) %>% Reduce(`+`, .)
      afull <- matrix(a10 - a20, nrow = p)

      if(length(ix) < p){
        R <- Rfull[ix,ix]
        R12 <- Rfull[ix,-ix]
        a <- afull[ix] - R12 %*% t(dat$f$fbar[j,-ix,drop = FALSE])
      }else{
        R <- Rfull
        a <- afull
      }
      S <- solve(R + T0)
      mu <- S %*% a
      return(list(Rfull = R, a = a, S = S, mu = mu))
    })
  }else{
    Oj <- map(dat$omega, function(o){o[j,]}) %>% unlist() %>%
      matrix(nrow = n, byrow = TRUE)
    samps <- lapply(seq(ns), function(i){
      asamp <- map(asamps, function(a){a[i,]}) %>% do.call(cbind, .)
      lsamp <- asamp %*% t(dat$G)
      Astar <- lapply(seq(p), function(jj){ # this is a list of W^{(a,j)}
        A <- t(asamp * Oj[,jj]) %*% asamp
        dat$G %*% A %*% t(dat$G)
      })
      Rfull <- Astar[[j]]
      a10 <- colSums(lsamp *rowSums(dat$Y*Oj))
      a20 <- lapply(seq(p)[-j], function(jj){
        Astar[[jj]]%*% t(dat$f$fbar[jj,,drop = FALSE])
      }) %>% Reduce(`+`, .)
      afull <- matrix(a10 - a20, nrow = p)

      if(length(ix) < p){
        R <- Rfull[ix,ix]
        R12 <- Rfull[ix,-ix]
        a <- afull[ix] - R12 %*% t(dat$f$fbar[j,-ix,drop = FALSE])
      }else{
        R <- Rfull
        a <- afull
      }
      S <- solve(R + T0)
      mu <- S %*% a
      return(list(R = R, a = a, S = S, mu = mu))
    })
  }
  mu_mat <- map(samps, "mu") %>% Reduce(cbind, .)

  mu_mean <- (map(samps, "mu") %>% Reduce(`+`, .))/ns
  Smean <- (map(samps, "S") %>% Reduce(`+`, .))/ns

  mu_mp <- mu_mean
  v_mp <- Smean + var(t(mu_mat))
  return(list(m = mu_mp, S = v_mp) )
}


update_beta_full_joint_mp <- function(dat, asamps){

  p <- dat$p
  n <- dat$n
  k <- dat$k
  ns <- nrow(asamps[[1]])

  ix <- p*(dat$beta$beta_k-1) + dat$beta$beta_j
  ix <- ix[!dat$beta$fix_beta]
  m <- length(ix)

  prior_precision <- dat$beta$prior_precision
  if(is.null(prior_precision)){
    T0 <- matrix(0, nrow = m, ncol = m)
  }else{
    T0 <- prior_precision
  }

  if(dat$s_equal){
    OYt <- dat$omega %*% t(dat$Y)
    samps <- lapplay(seq(ns), function(r){
      asamp <- map(asamps, function(a){a[r,]}) %>% do.call(cbind, .)
      lsamp <- asamp %*% t(dat$G)
      A <- t(asamp) %*% asamp
      Astar <- dat$G %*% A %*% t(dat$G)
      Rfull <- kronecker(Astar, dat$omega)

      afull <- lapply(seq(n), function(i){
        kronecker( matrix(lsamp[i,], nrow = k), matrix(OYt[,i], nrow = p))
      }) %>% Reduce(`+`, .)
      if(length(ix) < p*k){
        R <- Rfull[ix,ix]
        R12 <- Rfull[ix,-ix]
        fb <- matrix(as.vector(dat$f$fbar), ncol = 1)
        a <- afull[ix] - R12 %*% fb[-ix,,drop = F]
      }else{
        R <- Rfull
        a <- afull
      }
      S <- solve(R + T0)
      mu <- S %*% a
      return(list(S = S, mu = mu))
    })
  }else{
    OY <- lapply(seq(n), function(i){
      dat$omega[[i]] %*% matrix(dat$Y[i,], nrow = p)
    })
    samps <- lapply(seq(ns), function(r){
      asamp <- map(asamps, function(a){a[r,]}) %>% do.call(cbind, .)
      lsamp <- asamp %*% t(dat$G)

      Rfull <- lapply(seq(n), function(i){
        l <- asamp[i,]
        a <- outer(l, l)
        kronecker(tcrossprod(dat$G, tcrossprod(dat$G, a)), dat$omega[[i]]) ## kronecker(G %*% a %*% t(G), O)
      }) %>% Reduce(`+`, .)
      afull <- lapply(seq(n), function(i){
        kronecker(lsamp[i,], OY[[i]])
      }) %>% Reduce(`+`, .)
      if(length(ix) < p*k){
        R <- Rfull[ix,ix]
        R12 <- Rfull[ix,-ix]
        fb <- matrix(as.vector(dat$f$fbar), ncol = 1)
        a <- afull[ix] - R12 %*% fb[-ix,,drop = F]
      }else{
        R <- Rfull
        a <- afull
      }
      S <- solve(R + T0)
      mu <- S %*% a
      return(list(S = S, mu = mu))
    })
  }
  mu_mat <- map(samps, "mu") %>% Reduce(cbind, .)

  mu_mean <- (map(samps, "mu") %>% Reduce(`+`, .))/ns
  Smean <- (map(samps, "S") %>% Reduce(`+`, .))/ns

  mu_mp <- mu_mean
  v_mp <- Smean + var(t(mu_mat))
  return(list(m = mu_mp, S = v_mp) )
}

