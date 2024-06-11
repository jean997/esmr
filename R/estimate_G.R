#'@export
estimate_G <- function(beta_hat_X, se_X, R=NULL ,
                       type = c("gfa", "svd"),
                       svd_zthresh = 0,
                       augment = FALSE,
                       single_trait_thresh = 0.95,
                       add_trait1 = TRUE){
  type <- match.arg(type, type)
  beta_hat_X <- check_matrix(beta_hat_X)
  n <- nrow(beta_hat_X)
  p <- ncol(beta_hat_X)
   se_X <- check_matrix(se_X, n, p)
  #Z <- beta_hat_X/se_X
  #trait_scale <- c(1, apply(se_X[, -1, drop = FALSE]/se_X[,1], 2, median))
  colnames(beta_hat_X) <- colnames(se_X) <-  NULL

  if(is.null(R)){
    if(type == "gfa"){
      gfit <- GFA::gfa_fit(B_hat = beta_hat_X, S = se_X)
      myG <- gfit$F_hat
    }else if(type == "svd"){
      Z <- beta_hat_X/se_X
      Z[abs(Z) < svd_zthresh] <- 0
      myG <- svd(Z)$v
    }
  }else{
    R <- check_matrix(R, p, p)
    R <- check_R(R)
    if(type == "gfa"){
      gfit <- GFA::gfa_fit(B_hat = beta_hat_X, S = se_X, R = R)
      myG <- gfit$F_hat
    }else if(type == "svd"){
      Z <- beta_hat_X/se_X
      Z[abs(Z) < svd_zthresh] <- 0
      eR <- eigen(R)
      UTZ <-   Z %*% diag(1/sqrt(eR$values)) %*% eR$vectors
      V <- svd(UTZ)$v
      UV <- eR$vectors %*% diag(sqrt(eR$values)) %*% V
      myG <- UV
    }
  }
  if(augment){
    # Second augment version
    if(ncol(myG) < p){
      n_add <- p-ncol(myG)
      r <- rowSums(myG^2)
      which_add <- order(r)[1:n_add]
      A <- matrix(0, nrow = p, ncol = n_add)
      for(j in seq_along(which_add)) A[which_add[j],j] <- 1
      myG <- cbind(myG, A)
    }
  }
  #myG <- myG*trait_scale
  #myG <- sumstatFactors:::norm_cols(myG)$A
  k <- ncol(myG)
  if(!add_trait1) return(myG)
  G <- rbind(c(1, rep(0, k)),
             cbind(rep(0, p), myG))
  return(G)
}
