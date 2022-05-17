#'@export
var_beta <- function(Y, beta_m, lbar, l2bar, omega){
  n <- nrow(Y)
  p <- ncol(Y)
  stopifnot(nrow(lbar) == n & ncol(lbar) == p)
  stopifnot(nrow(l2bar) == n & ncol(l2bar) == p)

  if("matrix" %in% class(omega)){
    s_equal <- TRUE
  }else{
    stopifnot(class(omega) != "list")
    stopifnot(length(omega) == n)
    s_equal <- FALSE
  }
  beta_m <- as.vector(beta_m)

  stopifnot(length(beta_m) == p-1)

  s2l <- l2bar - (lbar^2)
  if(s_equal){
    ## Meat calculation (\hat{B}_n(\hat{\theta}_n))
    G <- t(t(Y - lbar)*omega[1,])
    gi <- rowSums(G)
    # E_lij_gi <- gi*lbar[,-1, drop =FALSE] # n by p-1
    # E_lij_gi <- E_lij_gi - s2l[,-1]*omega[1,-1]
    #
    # Bn <- lapply(seq(n), function(i){
    #   LL <- (beta_m*t(lbar[i,-1,drop = FALSE])) %*% lbar[i,-1,drop = FALSE]
    #   diag(LL) <- beta_m*l2bar[i,-1]
    #   beta_ll <- colSums(LL)*omega[1,1]
    #   grad <- t(E_lij_gi[i,,drop = FALSE]) - beta_ll
    #   ggt <- grad %*% t(grad)
    # }) %>% reduce(., `+`)
    # Bn <- Bn/n

    Bn <- lapply(seq(n), function(i){
      grad <- -(lbar[i,-1]*omega[1,1])*sum(beta_m*lbar[i,-1]) - beta_m*s2l[i,-1]*omega[1,1] +
               lbar[i,-1]*gi[i] - s2l[i,-1]*omega[1,-1]
      grad <- matrix(grad, nrow = p-1)
      ggt <- grad %*% t(grad)
    }) %>% reduce(., `+`)
    Bn <- Bn/n

    ## Bread Calculation
    bread1 <- lapply(seq(n), function(i){
      LL <- (t(lbar[i,-1,drop = FALSE]) %*% lbar[i,-1,drop = FALSE]) + diag(s2l[i,-1], nrow = p-1)
      -1*LL*omega[1,1]
    })
    b_m <- matrix(beta_m, nrow= p-1)
    blbar <- t(t(lbar[,-1,drop=FALSE])*beta_m)
    blbari <- rowSums(blbar)
    D2mm <- -(b_m %*% t(b_m))*omega[1,1] - b_m %*% omega[1,-1,drop = FALSE] - t(omega[1,-1,drop=FALSE])%*% t(b_m) - omega[-1, -1]
    D2m1j <- - c(0, beta_m)*omega[1,1] -omega[1,]
    D2mm <- rbind(D2m1j, cbind(D2m1j[-1], D2mm))
    bread2 <- lapply(seq(n), function(i){
      D2betam <- -omega[1,1]*t(b_m %*% lbar[i,-1,drop = FALSE]) -
                        t(t(omega[1,-1,drop = FALSE]) %*% lbar[i,-1,drop = FALSE])
      D2betam <- D2betam + (diag(p-1)*(gi[i] - omega[1,1]*blbari[i]))
      D2betam <- cbind(-t(lbar[i,-1, drop = FALSE])*omega[1,1], D2betam)
      D2betas <- diag( -beta_m*omega[1,1] - omega[1,-1], nrow = p-1)
      D2betas <- cbind(matrix(0, nrow = p-1, ncol = 1), D2betas)
      D2betap <- cbind(D2betam, D2betas)
      D2ss <- diag(-0.5*(1/(s2l[i,]^2)))
      zero <- matrix(0, nrow = p, ncol =p)
      D2pp <- rbind(cbind(D2mm, zero), cbind(zero, D2ss))
      D2betap %*% solve(D2pp) %*% t(D2betap)
    })

    Hn <- lapply(seq(n), function(i){
      bread1[[i]] - bread2[[i]]
    })%>% reduce(., `+`)
    Hn <- Hn/n
    Hninv <- solve(Hn)
    sand_var <- (1/n)*(Hninv %*% Bn %*% Hninv)
  }
  return(sand_var)
}

normal_means_loglik_d1 <- function(sample, g){

}

normal_means_lik_d1 <- function(sample, g){
  ix <- which(sample != 0)
  n <- length(sample)
  sample <- sample[ix]
  log_components1 <- sapply(sample, function(x){
    dnorm(x, mean = 0, sd = g$sd, log = TRUE )
  }) %>% t()
  log_components1 <- t(t(log_components1) + log(g$pi))
  components2 <- sapply(sample, function(x){
    dnorm(x, mean = 0, sd = g$sd) *(-x/(g$sd^2))*g$pi
  }) %>% t()
  res_log_denom <- apply(log_components1, 1, matrixStats::logSumExp )
  res_num <- apply(components2[,-1, drop = FALSE],1, sum)
  res <- res_num/exp(res_log_denom)
  return(sum(res)/n)
}

normal_means_lik_d2 <- function(sample, g){
  ix <- which(sample != 0)
  n <- length(sample)
  sample <- sample[ix]
  log_components1 <- sapply(sample, function(x){
    dnorm(x, mean = 0, sd = g$sd, log = TRUE )
  }) %>% t()
  log_components1 <- t(t(log_components1) + log(g$pi))
  components2 <- sapply(sample, function(x){
    dnorm(x, mean = 0, sd = g$sd) *(-x/(g$sd^2))*g$pi
  }) %>% t()
  components3 <- sapply(sample, function(x){
    dnorm(x, mean = 0, sd = g$sd) *((-1/(g$sd^2)) + x^2/(g$sd^2))*g$pi
  }) %>% t()
  res_log_denom <- apply(log_components1, 1, matrixStats::logSumExp )
  res_num1 <- apply(components2[,-1, drop = FALSE],1, sum)^2
  res_num2 <- apply(components3[,-1, drop = FALSE],1, sum)
  res_denom <- exp(res_log_denom)
  res <- (res_num1 - res_denom*res_num2)/res_denom
  return(sum(res)/n)
}
