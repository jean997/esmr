## likelihood function, very slow
log_py <- function(fit, g_hat, fbar, max_prob = 1, nmax = Inf){
  if(missing(g_hat)){
    g_hat <- fit$l$g_hat
  }
  if(missing(fbar)){
    fbar <- fit$f$fbar
  }
  fgbar <-  fbar %*% fit$G

  lpi <- lapply(g_hat, function(x){log(x$pi)})
  s <- lapply(g_hat, function(x){x$sd})

  if(max_prob < 1 | nmax < Inf){
    message("Identifying likelihood components.\n")
    lpi_mat <- unlist(lpi) %>% matrix(nrow = fit$p, byrow = T)
    top_combs <- get_top_combinations(x = lpi_mat, max_logsumexp = log(max_prob), nmax = nmax)
    m <- nrow(top_combs$combs)

    mmax <- sapply(lpi, length) %>% Reduce(`*`, .)

    lpi <- top_combs$values
    total_prob <- sum(exp(lpi))

    s_mat <- unlist(s) %>% matrix(nrow = fit$p, byrow = T)
    V <- apply(top_combs$combs, 1, function(c){
      ss <- s_mat[cbind(1:fit$p, c)]
      return(crossprod(t(fgbar)*ss, t(fgbar)*ss))
      #fgbar %*% diag(ss^2) %*% t(fgbar)
    }, simplify = F)
  }else{
    LPi <- expand.grid(lpi)
    lpi <- apply(LPi, 1, sum)
    S <- expand.grid(s)
    V <- apply(S, 1, function(s){
      crossprod(t(fgbar)*s, t(fgbar)*s)
      #fgbar %*% diag(s^2) %*% t(fgbar)
    }, simplify = F)
    m <- mmax <- length(V)
    total_prob <- 1
  }
  message(paste0("I will calculate ", m, " out of ", mmax, " total components of the likelihood, total integral ", total_prob, ".\n"))


  #equal_omega <- check_equal_omega(fit$omega)
  a <- -(fit$p/2)*log(2*base::pi)
  if(fit$s_equal){
    somega <- solve(fit$omega)
    lprob <- lapply(seq(m), function(mm){
      myV <- V[[mm]] + somega
      smV <- solve(myV)
      c <- -0.5* as.numeric(determinant(myV, logarithm = TRUE)$modulus)

      d <- sapply(1:fit$n, function(i){
        crossprod(crossprod(smV, fit$Y[i,]), fit$Y[i,])
      })

      lpi[mm] + c -0.5*d
    })
    A <- matrix(unlist(lprob) + a, nrow = fit$n, byrow = F)
    lp <- apply(A, 1, matrixStats::logSumExp)
    return(sum(lp))
  }else{
    lp <- sapply(1:fit$n, function(i){
      somega <- solve(fit$omega[[i]])
      lprob <- sapply(seq(m), function(mm){
        myV <- V[[mm]] + somega
        smV <- solve(myV)
        c <- - 0.5* as.numeric(determinant(myV, logarithm = TRUE)$modulus)

        d <- crossprod(crossprod(smV, fit$Y[i,]), fit$Y[i,])

        lpi[mm] + c -0.5*d
      }) + a
      matrixStats::logSumExp(lprob)
    })
    return(sum(lp))
  }
}






## gradient of log p(y | params)
grad_log_py <- function(fit, fbar, ix = NULL, max_prob = 1, nmax = Inf){ #Y, ghat, G, fbar, omega){
  #n <- nrow(fit$Y)
  #k <- ncol(fit$G)
  #p <- ncol(fit$Y)
  if(missing(fbar)){
    fbar <- fit$f$fbar
  }
  fgbar <- fbar %*% fit$G

  if(!is.null(ix)) fit <- subset_data(fit, ix)

  lpi <- lapply(fit$l$g_hat, function(x){log(x$pi)})
  s <- lapply(fit$l$g_hat, function(x){x$sd})


  if(max_prob < 1 | nmax < Inf){
    message("Identifying likelihood components.\n")
    lpi_mat <- unlist(lpi) %>% matrix(nrow = fit$p, byrow = T)
    top_combs <- get_top_combinations(x = lpi_mat, max_logsumexp = log(max_prob), nmax = nmax)
    m <- nrow(top_combs$combs)
    mmax <- sapply(lpi, length) %>% Reduce(`*`, .)

    lpi <- top_combs$values
    total_prob <- sum(exp(lpi))

    s_mat <- unlist(s) %>% matrix(nrow = fit$p, byrow = T)
    V <- apply(top_combs$combs, 1, function(c){
      s <- s_mat[cbind(1:fit$p, c)]
      crossprod(t(fgbar)*s, t(fgbar)*s)
      #fgbar %*% diag(s^2) %*% t(fgbar)
    }, simplify = FALSE)
    B <- apply(top_combs$combs, 1, function(c){
      s <- s_mat[cbind(1:fit$p, c)]
      crossprod(t(fit$G)*s, t(fit$G)*s)
    }, simplify = F)
  }else{
    LPi <- expand.grid(lpi)
    lpi <- apply(LPi, 1, sum)
    S <- expand.grid(s)
    V <- apply(S, 1, function(s){
      crossprod(t(fgbar)*s, t(fgbar)*s)
      #fgbar %*% diag(s^2) %*% t(fgbar)
    }, simplify = F)
    B <- apply(S, 1, function(s){
      crossprod(t(fit$G)*s, t(fit$G)*s)
    }, simplify = F)
    m <- mmax <- length(V)
    total_prob <- 1
  }
  message(paste0("I will calculate ", m, " out of ", mmax, " total components of the likelihood, total integral ", total_prob, ".\n"))


  ## caclulate d/df_ij V(F) for each variance matrix

  b_j <- fit$beta$beta_j[!fit$beta$fix_beta]
  b_k <- fit$beta$beta_k[!fit$beta$fix_beta]
  nvars <- length(b_j)
  E <- matrix(0, nrow = fit$p, ncol = fit$k)
  dV <- lapply(seq(nvars), function(i){
    myE <- E
    myE[b_j[i], b_k[i]] <- 1
    lapply(seq(m), function(mm){
      myE %*% B[[mm]] %*% t(fbar) + fbar %*% B[[mm]] %*% t(myE)
    })
  })
  # dV_is_nonzero <- lapply(seq(nvars), function(i){
  #   sapply(seq(m), function(mm){
  #     !all(dV[[i]][[mm]] == 0)
  #   }) %>% which()
  # })

  a <- -(fit$p/2)*log(2*base::pi)
  total_elts <- nvars + 1
  if(fit$s_equal){
    somega <- solve(fit$omega)
    wmod1 <- array(0, dim = c(fit$n, nvars))
    lprob <- sapply(seq(m), function(mm){
      myV <- V[[mm]] + somega
      smV <- solve(myV)
      c <- a - 0.5* as.numeric(determinant(myV, logarithm = TRUE)$modulus)

      d <- sapply(1:fit$n, function(i){
        crossprod(crossprod(smV, fit$Y[i,]), fit$Y[i,])
      })
      e <- c -0.5*d ## this part is N(y_i; 0, Sigma_k(beta) + S_i)

      ## calculate new weights for numerator of derivative
      for(nv in seq(nvars)){
        # if(!mm %in% dV_is_nonzero[[nv]]){
        #   wmod1[,nv] <- 0
        # }else{
        x1 <- smV %*% dV[[nv]][[mm]]
        x2 <- x1 %*% smV
        wmod1[,nv] <- -0.5*sum(diag(x1)) + 0.5*colSums(crossprod(x2, t(fit$Y)) *t(fit$Y))
        #}
      }
      # wmod1 <- sapply(seq(nvars), function(nv){
      #   if(!mm %in% dV_is_nonzero[[nv]]) return(rep(0, fit$n))
      #   x1 <- smV %*% dV[[nv]][[mm]]
      #   x2 <- x1 %*% smV
      #   sapply(1:fit$n, function(i){-0.5*sum(diag(x1)) + 0.5*crossprod(crossprod(x2, fit$Y[i,]), fit$Y[i,])})
      # }, simplify = "array")
      #lpi[mm] + c -(1/2)*d
      return(cbind(e, wmod1))
    }, simplify = "array")

    #A <- matrix(unlist(lprob), nrow = n, byrow = F)
    #lp <- apply(A, 1, matrixStats::logSumExp)
    x <- t(t(lprob[,1,]) + lpi)
    lf <- apply(x, 1, matrixStats:::logSumExp)
    #fdot <- sum(lprob[2,]*exp(x))
    isneg <- lprob[,-1,] < 0
    lfdotpos <- apply(lprob[,-1,]*(!isneg), 2, function(xx){
      apply(log(xx) + x, 1, matrixStats:::logSumExp)
    })
    lfdotneg <- apply(lprob[,-1,]*(isneg), 2, function(xx){
      apply(log(-xx) + x, 1, matrixStats:::logSumExp)
    })
    fdot <- exp(lfdotpos) - exp(lfdotneg)
    #c(lf, fdot/exp(lf))

    gr <- sign(fdot)*exp(log(abs(fdot)) - lf) # fdot/exp(lf)
    In <- lapply(seq(fit$n), function(i){
      outer(gr[i,], gr[i,])
    }) %>% Reduce(`+`, .)
    In <- In/fit$n
    Sn <- colSums(gr)/fit$n
    return(list(log_py = sum(lf), grad = Sn*fit$n,
                Sn = Sn,
                In = In))
  }else{
    wmod1 <- numeric(nvars)
    lprob <- array(0, dim = c(total_elts, m))
    #lp <- array(0, dim = c(total_elts, fit$n))
    lp <- sapply(1:fit$n, function(i){
    #for(i in seq(fit$n)){
      somega <- solve(fit$omega[[i]])

      #lprob <- sapply(seq(m), function(mm){
      for(mm in seq(m)){
        myV <- V[[mm]] + somega
        smV <- solve(myV)
        c <- - 0.5* as.numeric(determinant(myV, logarithm = TRUE)$modulus)

        d <- crossprod(crossprod(smV, fit$Y[i,]), fit$Y[i,]) %>% as.numeric
        e <- c -0.5*d ## this part is N(y_i; 0, Sigma_k(beta) + S_i), will add a later

        ## calculate new weights for numerator of derivative
        #wmod1 <- sapply(seq(nvars), function(nv){

        for(nv in seq(nvars)){
          #if(!mm %in% dV_is_nonzero[[nv]]) return(0)
          # x1 <- smV %*% dV[[nv]][[mm]]
          # x2 <- x1 %*% smV
          # wmod1[nv] <-
          #     -0.5*sum(diag(x1)) + 0.5*crossprod(crossprod(x2, fit$Y[i,]), fit$Y[i,])

          t1 <- -0.5*sum(smV * t(dV[[nv]][[mm]])) # -0.5*Trace(smV %*% dV[[nv]][[mm]])
          smV_Y <- crossprod(smV, fit$Y[i,])
          t2 <- 0.5*crossprod(crossprod(dV[[nv]][[mm]], smV_Y), smV_Y)
          wmod1[nv] <- t1 + t2
        }#)

        #return(c(e, wmod1)) #, as.vector(wmod2)))
        lprob[,mm] <- c(e, wmod1)
      }# , simplify =  "array")
      #matrixStats::logSumExp(lprob)
      x <- lprob[1,] + lpi + a
      lf <- matrixStats:::logSumExp(x)  ## log(f(beta; x_i))

      #fdot <- sum(lprob[2,]*exp(x))
      isneg <- t(lprob[-1,,drop = F] < 0)
      lfdotpos <- apply(t(lprob[-1,,drop = F])*(!isneg), 2, function(xx){
        matrixStats::logSumExp(log(xx) + x)
      })
      lfdotneg <- apply(t(lprob[-1,,drop = F])*(isneg), 2,  function(xx){
        matrixStats::logSumExp(log(-1*xx) + x)
      })
      fdot <- exp(lfdotpos) - exp(lfdotneg) ## grad f(beta; x_i) and also second derivatives
      grad <- sign(fdot)*exp(log(abs(fdot)) - lf)#fdot/exp(lf)
      #lp[, mm] <-
      c(lf, grad)
    }, simplify = "array")
    X <- rowSums(lp)

    In <- lapply(seq(fit$n), function(i){
      outer(lp[-1, i], lp[-1, i])
    }) %>% Reduce(`+`, .)
    grad = X[1 + seq(nvars)]

    return(list(log_py = X[1], grad = grad,
                Sn = grad/fit$n,
                In = In/fit$n))
  }
}


hess_log_py <- function(fit, fbar, ix = NULL, max_prob = 1, nmax = Inf){

  if(missing(fbar)){
    fbar <- fit$f$fbar
  }


  fgbar <- fbar %*% fit$G

  if(!is.null(ix)) fit <- subset_data(fit, ix)

  lpi <- lapply(fit$l$g_hat, function(x){log(x$pi)})
  s <- lapply(fit$l$g_hat, function(x){x$sd})

  if(max_prob < 1 | nmax < Inf){
    message("Identifying likelihood components.\n")
    lpi_mat <- unlist(lpi) %>% matrix(nrow = fit$p, byrow = T)
    top_combs <- get_top_combinations(x = lpi_mat, max_logsumexp = log(max_prob), nmax = nmax)
    m <- nrow(top_combs$combs)
    mmax <- sapply(lpi, length) %>% Reduce(`*`, .)

    lpi <- top_combs$values
    total_prob <- sum(exp(lpi))

    s_mat <- unlist(s) %>% matrix(nrow = fit$p, byrow = T)
    V <- apply(top_combs$combs, 1, function(c){
      s <- s_mat[cbind(1:fit$p, c)]
      crossprod(t(fgbar)*s, t(fgbar)*s)
      #fgbar %*% diag(s^2) %*% t(fgbar)
    }, simplify = FALSE)
    B <- apply(top_combs$combs, 1, function(c){
      s <- s_mat[cbind(1:fit$p, c)]
      crossprod(t(fit$G)*s, t(fit$G)*s)
    }, simplify = F)
  }else{
    LPi <- expand.grid(lpi)
    lpi <- apply(LPi, 1, sum)
    S <- expand.grid(s)
    V <- apply(S, 1, function(s){
      crossprod(t(fgbar)*s, t(fgbar)*s)
      #fgbar %*% diag(s^2) %*% t(fgbar)
    }, simplify = F)
    B <- apply(S, 1, function(s){
      crossprod(t(fit$G)*s, t(fit$G)*s)
    }, simplify = F)
    m <- mmax <- length(V)
    total_prob <- 1
  }

  ## caclulate d/df_ij V(F) for each variance matrix
  nvars <- length(fit$beta$beta_j)
  E <- matrix(0, nrow = fit$p, ncol = fit$k)
  dV <- lapply(seq(nvars), function(i){
    myE <- E
    myE[fit$beta$beta_j[i], fit$beta$beta_k[i]] <- 1
    lapply(seq(m), function(mm){
      myE %*% B[[mm]] %*% t(fbar) + fbar %*% B[[mm]] %*% t(myE)
    })
  })
  dV_is_nonzero <- lapply(seq(nvars), function(i){
    sapply(seq(m), function(mm){
      !all(dV[[i]][[mm]] == 0)
    }) %>% which()
  })

  d2V <- lapply(seq(nvars), function(i){
    myEi <- E
    myEi[fit$beta$beta_j[i], fit$beta$beta_k[i]] <- 1
    lapply(seq(nvars), function(j){
      myEj <- E
      myEj[fit$beta$beta_j[j], fit$beta$beta_k[j]] <- 1
      lapply(seq(m), function(mm){
        myEi %*% B[[mm]] %*% t(myEj) + myEj %*% B[[mm]] %*% t(myEi)
      })
    })
  })

  #equal_omega <- check_equal_omega(fit$omega)
  a <- -(fit$p/2)*log(2*base::pi)
  total_elts <- 1 + nvars + nvars^2
  if(fit$s_equal){
    somega <- solve(fit$omega)
    wmod1 <- array(0, dim = c(fit$n, nvars))
    wmod2 <- array(0, dim = c(fit$n, nvars, nvars))
    lprob <- array(0, dim = c(fit$n, total_elts, m))
    #lprob <- sapply(seq(m), function(mm){
    for(mm in seq(m)){
      myV <- V[[mm]] + somega
      smV <- solve(myV)
      c <- a - 0.5* as.numeric(determinant(myV, logarithm = TRUE)$modulus)
      d <- colSums(crossprod(smV, t(fit$Y)) * t(fit$Y))
      # d <- sapply(1:fit$n, function(i){
      #   crossprod(crossprod(smV, fit$Y[i,]), fit$Y[i,])
      # })
      e <- c -0.5*d ## log N(y_i; 0, myV)

      smV_dV <- lapply(seq(nvars), function(nv){
        smV %*% dV[[nv]][[mm]]
      })

      ## calculate new weights for numerator of derivative
      #wmod1 <- sapply(seq(nvars), function(nv){
      for(nv in seq(nvars)){
        if(!mm %in% dV_is_nonzero[[nv]]){
          wmod1[,nv] <- 0
        }else{
        #x1 <- smV %*% dV[[nv]][[mm]]
          x2 <- smV_dV[[nv]] %*% smV
          wmod1[,nv] <- -0.5*sum(diag(smV_dV[[nv]])) + 0.5*colSums(crossprod(x2, t(fit$Y)) *t(fit$Y))
        }
      }#, simplify = "array")

      #wmod2 <- sapply(seq(nvars), function(nv1){
        #x_nv1 <- smV %*% dV[[nv1]][[mm]]
      #  sapply(seq(nvars), function(nv2){
      for(nv1 in seq(nvars)){
        for(nv2 in seq(nvars)){
          #if(!mm %in% d2V_is_nonzero[[nv1]][[nv2]]) return(0)
          #x_nv2 <- smV %*% dV[[nv2]][[mm]]

          x1 <- dV[[nv2]][[mm]]%*%smV_dV[[nv1]] - d2V[[nv1]][[nv2]][[mm]]
          x2 <- -smV %*% x1
          x3 <- x1 + dV[[nv1]][[mm]]%*%smV_dV[[nv2]]
          x4 <- -smV %*% x3 %*% smV

          wmod2[,nv1, nv2] <- -0.5*sum(diag(x2)) +0.5*colSums(crossprod(x4, t(fit$Y)) *t(fit$Y))
        }
      }#, simplify =  "array")
      for(i in seq(fit$n)){
        wmod2[i,,] <- wmod2[i,,] + tcrossprod(wmod1[i,])
      }
      wmod2_long <- matrix(wmod2, ncol = nvars^2, nrow = fit$n )
      lprob[,,mm] <- cbind(e, wmod1, wmod2_long)
    }#, simplify = "array")

    x <- t(t(lprob[,1,]) + lpi)
    lf <- apply(x, 1, matrixStats:::logSumExp)

    isneg <- lprob[,-1,] < 0
    lfdotpos <- apply(lprob[,-1,]*(!isneg), 2, function(xx){
      apply(log(xx) + x, 1, matrixStats:::logSumExp)
    })
    lfdotneg <- apply(lprob[,-1,]*(isneg), 2, function(xx){
      apply(log(-xx) + x, 1, matrixStats:::logSumExp)
    })

    fdots <- exp(lfdotpos) - exp(lfdotneg) ## grad f(beta; x_i) and also second derivatives
    grad <- colSums(fdots[,seq(nvars)]/exp(lf))

    hess <- lapply(1:fit$n, function(i){
      tcrossprod(fdots[i,seq(nvars)])/(exp(lf[i])^2) - matrix(fdots[i, -seq(nvars)], nrow = nvars)/exp(lf[i])
    }) %>% Reduce("+", .)

    return(list(log_py = sum(lf), grad = grad,
                hess = hess))
  }else{
    wmod1 <- numeric(nvars)
    wmod2 <- array(0, dim = c(nvars, nvars))
    lprob <- array(0, dim = c(total_elts, m))
    #lp <- array(0, dim = c(total_elts, fit$n))
    lp <- sapply(seq(fit$n), function(i){
    #for(i in seq(fit$n)){
      somega <- solve(fit$omega[[i]])
      #lprob <- sapply(seq(m), function(mm){
      for(mm in seq(m)){
        myV <- V[[mm]] + somega
        smV <- solve(myV)
        c <- a - 0.5* as.numeric(determinant(myV, logarithm = TRUE)$modulus)

        d <- crossprod(crossprod(smV, fit$Y[i,]), fit$Y[i,]) %>% as.numeric
        e <- c -0.5*d ## this part is N(y_i; 0, Sigma_k(beta) + S_i)

        smV_dV <- lapply(seq(nvars), function(nv){
          smV %*% dV[[nv]][[mm]]
        })

        ## calculate new weights for numerator of derivative
        # wmod1
        for (nv in seq(nvars)) {
          # if (!mm %in% dV_is_nonzero[[nv]]){
          #   wmod1[nv] <- 0
          # }
          x2 <- smV_dV[[nv]] %*% smV
          wmod1[nv] <- -0.5 * sum(diag(smV_dV[[nv]])) + 0.5 * as.numeric(crossprod(crossprod(x2, fit$Y[i,]), fit$Y[i,]))
        }
        # wmod2
        for (nv1 in seq(nvars)) {
          for (nv2 in seq(nvars)) {
            x1 <- dV[[nv2]][[mm]]%*%smV_dV[[nv1]] - d2V[[nv1]][[nv2]][[mm]]
            x2 <- -smV %*% x1
            x3 <- x1 + dV[[nv1]][[mm]]%*%smV_dV[[nv2]]
            x4 <- -smV %*% x3 %*% smV
            wmod2[nv1, nv2] <- -0.5*sum(diag(x2)) +0.5*crossprod(crossprod(x4, fit$Y[i,]), fit$Y[i,])
          }
        }
        wmod2 <- wmod2 + tcrossprod(wmod1)
        #a2 <- wmod1 + e
        #return(c(e, wmod1, as.vector(wmod2)))
        lprob[,mm] <- c(e, wmod1, as.vector(wmod2))
        #return(a1)
      }
      x <- lprob[1,] + lpi
      lf <- matrixStats:::logSumExp(x)  ## log(f(beta; x_i))

      isneg <- t(lprob[-1,] < 0)
      lfdotpos <- apply(t(lprob[-1,])*(!isneg), 2, function(xx){
        matrixStats::logSumExp(log(xx) + x)
      })
      lfdotneg <- apply(t(lprob[-1,])*(isneg), 2,  function(xx){
        matrixStats::logSumExp(log(-1*xx) + x)
      })
      fdots <- exp(lfdotpos) - exp(lfdotneg) ## grad f(beta; x_i) and also second derivatives
      fdot <- fdots[seq(nvars)]
      fdotdot <- matrix(fdots[-seq(nvars)], nrow = nvars)

      hess <- tcrossprod(fdot)/(exp(lf)^2) - fdotdot/exp(lf)
      grad <- fdot/exp(lf)
      c(lf, grad, as.vector(hess))
      #lp[,i] <- c(lf, grad, as.vector(hess))
    }, simplify = "array")
    X <- rowSums(lp)

    return(list(log_py = X[1], grad = X[1 + seq(nvars)],
                hess = matrix(X[-seq(nvars + 1)], nrow = nvars)))
  }
}


#'@title One step update
#'@param fit ESMR fit object
#'@param max_steps Number of steps to take (default is 1)
#'@param tol Tolerance for convergence if max_steps is > 1
#'@param calc_hess Calculate the hesssian? (this can be slow)
#'@param max_components Maximum number of likelihood components to use in approximation
#'@param max_prob Maximum integral of approximate likelihood
#'@param sub_size Size of subset to use to calculate the gradient of the likelihood.
#'@details This function computes a one step update from the variational solution provided by ESMR. By default
#'the gradient of the likelihood is computed exactly. This should work well if the number of traits is < 15.
#'For larger numbers of traits, you may want to use an approximation to the likelihood controlled by parameters max_components,
#'max_prob and sub_size.
#'@export
optimize_lpy2 <- function(fit,
                         max_steps = 1,
                         tol = 1e-5,
                         calc_hess = FALSE,
                         max_components = Inf,
                         max_prob = 1,
                         sub_size = fit$n){

  i <- 1
  bj <- fit$beta$beta_j
  bk <- fit$beta$beta_k
  update_fbar <- function(fbar, new_beta){
    for(ii in seq(bj)){
      fbar[bj[ii], bk[ii]] <- new_beta[ii]
    }
    return(fbar)
  }
  fbar <- fit$f$fbar
  beta <- fit$beta$beta_m
  done <- FALSE
  if(sub_size >= fit$n){
    ix <- NULL
  }
  while(i <= max_steps & !done){
    cat(i, " ")
    if(sub_size < fit$n){
      ix <- sample(seq(fit$n), size = sub_size, replace = FALSE)
      ix <- sort(ix)
    }
    g <- grad_log_py(fit, fbar, ix = ix,
                     max_prob = max_prob,
                     nmax = max_components)
    step <- solve(g$In) %*% g$Sn
    cat(step, " ")
    beta <- beta + step
    fbar <- update_fbar(fbar, beta)
    if(all(abs(step) < tol)) done <- TRUE
    i <- i + 1
    cat(g$log_py, "\n")
  }
  fit$beta$beta_m <- beta
  fit$f$fbar <- fbar
  fit$f$fgbar <- fit$G %*% fbar
  if(calc_hess){
    h <- hess_log_py(fit, fbar, ix = ix,
                     max_prob = max_prob,
                     nmax = max_components)
    fit$beta$V <- solve(h$hess)
    fit$beta$beta_s <- sqrt(diag(fit$beta$V))
    fit$direct_effects <- total_to_direct(t(fit$f$fbar) - diag(fit$p))
    delt_pvals <- delta_method_pvals(fit)
    fit$pvals_dm <- delt_pvals$pmat
    fit$se_dm <- delt_pvals$semat
    fit$likelihood <- h$log_py
  }else{
    fit$beta$V <- solve(g$In)/fit$n
    fit$beta$beta_s <- sqrt(diag(fit$beta$V))
    fit$direct_effects <- total_to_direct(t(fit$f$fbar) - diag(fit$p))
    delt_pvals <- delta_method_pvals(fit)
    fit$pvals_dm <- delt_pvals$pmat
    fit$se_dm <- delt_pvals$semat
    #fit$likelihood <- log_py(fit)
  }
  return(fit)
}




