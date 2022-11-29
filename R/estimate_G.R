#'@export
estimate_G <- function(beta_hat_X, se_X, R=NULL ,
                       type = c("gfa", "svd"), svd_zthresh = 0){
  type <- match.arg(type, type)
  beta_hat_X <- check_matrix(beta_hat_X, "beta_hat_X")
  n <- nrow(beta_hat_X)
  p <- ncol(beta_hat_X)
  se_X <- check_matrix(se_X, "se_X", n, p)
  Z <- beta_hat_X/se_X

  if(is.null(R)){
    if(type == "gfa"){
      gfit <- sumstatFactors::gfa_fit(Z_hat = Z)
      myG <- gfit$F_hat_scaled
    }else if(type == "svd"){
      Z[abs(Z) < svd_zthresh] <- 0
      myG <- svd(Z)$v
    }
  }else{
    R <- check_matrix(R, "R", p, p)
    R <- check_R(R)
    if(type == "gfa"){
      gfit <- sumstatFactors::gfa_fit(Z_hat = Z, R = R)
      myG <- gfit$F_hat_scaled
    }else if(type == "svd"){
      Z[abs(Z) < svd_zthresh] <- 0
      eR <- eigen(R)
      UTZ <-   Z %*% diag(1/sqrt(eR$values)) %*% eR$vectors
      V <- svd(UTZ)$v
      UV <- eR$vectors %*% diag(sqrt(eR$values)) %*% V
      myG <- UV
    }
  }
  k <- ncol(myG)
  G <- rbind(c(1, rep(0, k)),
             cbind(rep(0, p), myG))
  return(G)
}

pad_G <- function(Gx, p, xindex){
  Gx <- check_matrix(Gx, "Gx", length(xindex), length(xindex))
  G <- diag(p)
  G[xindex, xindex] <- Gx
  return(G)
}
