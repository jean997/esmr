#'@export
estimate_G <- function(beta_hat_X, se_X, R=NULL ,
                       type = c("gfa", "svd"),
                       svd_zthresh = 0,
                       augment = FALSE,
                       single_trait_thresh = 0.95){
  type <- match.arg(type, type)
  beta_hat_X <- check_matrix(beta_hat_X, "beta_hat_X")
  n <- nrow(beta_hat_X)
  p <- ncol(beta_hat_X)
  se_X <- check_matrix(se_X, "se_X", n, p)
  Z <- beta_hat_X/se_X

  trait_scale <- c(1, apply(se_X[, -1, drop = FALSE]/se_X[,1], 2, median))

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
  if(augment){
    # cmax <- apply(abs(myG), 2, max)
    # if(any(cmax > single_trait_thresh)){
    #   st_factors <- which(cmax > single_trait_thresh)
    #   which_st <- apply(abs(myG[,st_factors]), 2, which.max)
    #   if(!all((1:p) %in% which_st)){
    #     n_add <- sum(!(1:p) %in% which_st)
    #     which_add <- (1:p)[!(1:p) %in% which_st]
    #     A <- matrix(0, nrow = p, ncol = n_add)
    #     for(j in seq_along(which_add)) A[which_add[j],j] <- 1
    #     myG <- cbind(myG, A)
    #   }
    # }else{
    #   myG <- cbind(myG, diag(p))
    # }

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
  myG <- myG*trait_scale
  myG <- sumstatFactors:::norm_cols(myG)
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
