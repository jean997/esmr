## gradient of log p(y | params)
grad_log_py <- function(fit, fbar){ #Y, ghat, G, fbar, omega){
  #n <- nrow(fit$Y)
  #k <- ncol(fit$G)
  #p <- ncol(fit$Y)
  fgbar <- fbar %*% fit$G

  lpi <- lapply(fit$l$g_hat, function(x){log(x$pi)})
  LPi <- expand.grid(lpi)
  lpi <- apply(LPi, 1, sum)
  pi <- exp(lpi)

  s <- lapply(fit$l$g_hat, function(x){x$sd})
  S <- expand.grid(s)
  V <- apply(S, 1, function(s){
    fgbar %*% diag(s^2) %*% t(fgbar)
  }, simplify = F)

  m <- length(V)

  B <- apply(S, 1, function(s){
    fit$G %*% diag(s^2) %*% t(fit$G)
  }, simplify = F)

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

  #equal_omega <- check_equal_omega(fit$omega)

  if(fit$s_equal){
    somega <- solve(fit$omega)
    lprob <- sapply(seq(m), function(mm){
      myV <- V[[mm]] + somega
      smV <- solve(myV)
      c <- -(fit$p/2)*log(2*base::pi) - (1/2)* as.numeric(determinant(myV, logarithm = TRUE)$modulus)

      d <- sapply(1:fit$n, function(i){
        crossprod(crossprod(smV, fit$Y[i,]), fit$Y[i,])
      })
      e <- c -(1/2)*d ## this part is N(y_i; 0, Sigma_k(beta) + S_i)

      ## calculate new weights for numerator of derivative
      wmod1 <- sapply(seq(nvars), function(nv){
        if(!mm %in% dV_is_nonzero[[nv]]) return(rep(0, fit$n))
        x1 <- smV %*% dV[[nv]][[mm]]
        x2 <- x1 %*% smV
        sapply(1:fit$n, function(i){-0.5*sum(diag(x1)) + 0.5*crossprod(crossprod(x2, fit$Y[i,]), fit$Y[i,])})
      }, simplify = "array")
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

    return(list(log_py = sum(lf), grad = colSums(fdot/exp(lf))))
  }else{

    lp <- sapply(1:fit$n, function(i){
      somega <- solve(fit$omega[[i]])

      lprob <- sapply(seq(m), function(mm){
        myV <- V[[mm]] + somega
        smV <- solve(myV)
        c <- -(fit$p/2)*log(2*base::pi) - (1/2)* as.numeric(determinant(myV, logarithm = TRUE)$modulus)

        d <- crossprod(crossprod(smV, fit$Y[i,]), fit$Y[i,]) %>% as.numeric
        e <- c -(1/2)*d ## this part is N(y_i; 0, Sigma_k(beta) + S_i)

        ## calculate new weights for numerator of derivative
        wmod1 <- sapply(seq(nvars), function(nv){
          if(!mm %in% dV_is_nonzero[[nv]]) return(0)
          x1 <- smV %*% dV[[nv]][[mm]]
          x2 <- x1 %*% smV
          -0.5*sum(diag(x1)) + 0.5*crossprod(crossprod(x2, fit$Y[i,]), fit$Y[i,])
        })

        return(c(e, wmod1)) #, as.vector(wmod2)))
        #return(a1)
      }, simplify =  "array")
      #matrixStats::logSumExp(lprob)
      x <- lprob[1,] + lpi
      lf <- matrixStats:::logSumExp(x)  ## log(f(beta; x_i))

      #fdot <- sum(lprob[2,]*exp(x))
      isneg <- t(lprob[-1,] < 0)
      lfdotpos <- apply(t(lprob[-1,])*(!isneg), 2, function(xx){
        matrixStats::logSumExp(log(xx) + x)
      })
      lfdotneg <- apply(t(lprob[-1,])*(isneg), 2,  function(xx){
        matrixStats::logSumExp(log(-1*xx) + x)
      })
      fdot <- exp(lfdotpos) - exp(lfdotneg) ## grad f(beta; x_i) and also second derivatives
      grad <- fdot/exp(lf)
      c(lf, grad)
    }, simplify = "array")
    X <- rowSums(lp)

    return(list(log_py = X[1], grad = X[1 + seq(nvars)]))
  }
}


hess_log_py <- function(fit, fbar){
  fgbar <- fbar %*% fit$G

  lpi <- lapply(fit$l$g_hat, function(x){log(x$pi)})
  LPi <- expand.grid(lpi)
  lpi <- apply(LPi, 1, sum)
  pi <- exp(lpi)

  s <- lapply(fit$l$g_hat, function(x){x$sd})
  S <- expand.grid(s)
  V <- apply(S, 1, function(s){
    fgbar %*% diag(s^2) %*% t(fgbar)
  }, simplify = F)

  m <- length(V)

  B <- apply(S, 1, function(s){
    fit$G %*% diag(s^2) %*% t(fit$G)
  }, simplify = F)

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

  if(fit$s_equal){
    somega <- solve(fit$omega)
    lprob <- sapply(seq(m), function(mm){
      myV <- V[[mm]] + somega
      smV <- solve(myV)
      c <- -(fit$p/2)*log(2*base::pi) - (1/2)* as.numeric(determinant(myV, logarithm = TRUE)$modulus)

      d <- sapply(1:fit$n, function(i){
        crossprod(crossprod(smV, fit$Y[i,]), fit$Y[i,])
      })
      e <- c -(1/2)*d ## log N(y_i; 0, myV)

      ## calculate new weights for numerator of derivative
      wmod1 <- sapply(seq(nvars), function(nv){
        if(!mm %in% dV_is_nonzero[[nv]]) return(rep(0, fit$n))
        x1 <- smV %*% dV[[nv]][[mm]]
        x2 <- x1 %*% smV
        sapply(1:fit$n, function(i){-0.5*sum(diag(x1)) + 0.5*crossprod(crossprod(x2, fit$Y[i,]), fit$Y[i,])})
      }, simplify = "array")

      wmod2 <- sapply(seq(nvars), function(nv1){
        x_nv1 <- smV %*% dV[[nv1]][[mm]]
        sapply(seq(nvars), function(nv2){
          #if(!mm %in% d2V_is_nonzero[[nv1]][[nv2]]) return(0)
          x_nv2 <- smV %*% dV[[nv2]][[mm]]

          x1 <- dV[[nv2]][[mm]]%*%x_nv1 - d2V[[nv1]][[nv2]][[mm]]
          x2 <- -smV %*% x1
          x3 <- x1 + dV[[nv1]][[mm]]%*%x_nv2
          x4 <- -smV %*% x3 %*% smV

          sapply(1:fit$n, function(i){-0.5*sum(diag(x2)) +0.5*crossprod(crossprod(x4, fit$Y[i,]), fit$Y[i,])})
        })
      }, simplify =  "array")
      wmod2 <- sapply(1:fit$n, function(i){
        wmod2[i,,] + tcrossprod(wmod1[i,])
      }, simplify = "array")
      wmod2 <- t(matrix(wmod2, nrow = nvars^2, ncol = 5))

      return(cbind(e, wmod1, wmod2))
    }, simplify = "array")

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

    lp <- sapply(1:fit$n, function(i){
      somega <- solve(fit$omega[[i]])

      lprob <- sapply(seq(m), function(mm){
        myV <- V[[mm]] + somega
        smV <- solve(myV)
        c <- -(fit$p/2)*log(2*base::pi) - (1/2)* as.numeric(determinant(myV, logarithm = TRUE)$modulus)

        d <- crossprod(crossprod(smV, fit$Y[i,]), fit$Y[i,]) %>% as.numeric
        e <- c -(1/2)*d ## this part is N(y_i; 0, Sigma_k(beta) + S_i)

        ## calculate new weights for numerator of derivative
        wmod1 <- sapply(seq(nvars), function(nv){
          if(!mm %in% dV_is_nonzero[[nv]]) return(0)
          x1 <- smV %*% dV[[nv]][[mm]]
          x2 <- x1 %*% smV
          -0.5*sum(diag(x1)) + 0.5*crossprod(crossprod(x2, fit$Y[i,]), fit$Y[i,])
        })
        wmod2 <- sapply(seq(nvars), function(nv1){
          x_nv1 <- smV %*% dV[[nv1]][[mm]]
          sapply(seq(nvars), function(nv2){
            #if(!mm %in% d2V_is_nonzero[[nv1]][[nv2]]) return(0)
            x_nv2 <- smV %*% dV[[nv2]][[mm]]

            x1 <- dV[[nv2]][[mm]]%*%x_nv1 - d2V[[nv1]][[nv2]][[mm]]
            x2 <- -smV %*% x1
            x3 <- x1 + dV[[nv1]][[mm]]%*%x_nv2
            x4 <- -smV %*% x3 %*% smV

            -0.5*sum(diag(x2)) +0.5*crossprod(crossprod(x4, fit$Y[i,]), fit$Y[i,])
          })
        }, simplify =  "array")
        wmod2 <- wmod2 + tcrossprod(wmod1)
        #a2 <- wmod1 + e
        return(c(e, wmod1, as.vector(wmod2)))
        #return(a1)
      }, simplify =  "array")
      #matrixStats::logSumExp(lprob)
      x <- lprob[1,] + lpi
      lf <- matrixStats:::logSumExp(x)  ## log(f(beta; x_i))

      #fdot <- sum(lprob[2,]*exp(x))
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
    }, simplify = "array")
    X <- rowSums(lp)

    return(list(log_py = X[1], grad = X[1 + seq(nvars)],
                hess = matrix(X[-seq(nvars + 1)], nrow = nvars)))
  }
}

optimize_lpy <- function(fit,
                         max_steps = 10,
                         tol = 1e-5,
                         calc_hess = FALSE){
  V <- fit$beta$V
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
  while(i <= max_steps & !done){
    cat(i, " ")
    g <- grad_log_py(fit, fbar)
    step <- V %*% matrix(g$grad, ncol = 1)
    cat(step, "\n")
    beta <- beta + step
    fbar <- update_fbar(fbar, beta)
    if(all(abs(step) < tol)) done <- TRUE
    i <- i + 1
  }
  fit$beta$beta_m <- beta
  fit$f$fbar <- fbar
  fit$f$fgbar <- fit$G %*% fbar
  if(calc_hess){
    h <- hess_log_py(fit, fbar)
    fit$beta$V <- solve(h$hess)
    fit$beta$beta_s <- sqrt(diag(fit$beta$V))
    fit$direct_effects <- total_to_direct(t(fit$f$fbar) - diag(fit$p))
    delt_pvals <- delta_method_pvals(fit)
    fit$pvals_dm <- delt_pvals$pmat
    fit$se_dm <- delt_pvals$semat
    fit$likelihood <- h$log_py
  }
  return(fit)
}
