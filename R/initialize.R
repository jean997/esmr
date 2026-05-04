# previously init_l_future
init_l <- function(n, p, m){
  lbar <- matrix(0, nrow = n, ncol = p)
  abar <- matrix(0, nrow = n, ncol = m)
  lfsr <- matrix(1, nrow = n, ncol = m)
  g_hat <- list()
  return(list(lbar = lbar, l2bar = lbar,
              abar = abar, a2bar = abar,
              #wpost = lbar, mupost = lbar, s2post = lbar,
              #post_mode = abar,
              lfsr = lfsr, g_hat = g_hat))
}



init_beta <- function(dat, init_beta_X_Y = NULL){

  dat$beta <- list()

  if(dat$is_nesmr){
    B <- check_B_template(dat$B_template, dat$p, restrict_dag = dat$restrict_dag)
    which_beta <- rbind(B$which_tot_u, B$which_tot_c)[,c(2,1), drop=FALSE] ## transpose
    dat$beta$fix_beta <- c(rep(FALSE, nrow(B$which_tot_u)), rep(TRUE, nrow(B$which_tot_c)))
  }else if(dat$is_factors){
    which_beta <- cbind(c(rep(1,  dat$k - 1),
                          rep(2, dat$k-2)), c(2:dat$k , 3:dat$k))
    dat$beta$fix_beta <- rep(FALSE, 2*dat$k-3)
  }else{
    which_beta <- cbind(rep(1,  dat$p - 1), 2:dat$p )
    dat$beta$fix_beta <- rep(FALSE, dat$p-1)
  }

  colnames(which_beta) <- c("row", "col")
  dat$beta$beta_j <- which_beta[,1]
  dat$beta$beta_k <- which_beta[,2]

  nb <- length(dat$beta$beta_j)
  if(dat$is_nesmr){
    dat$beta$beta_m <- dat$B_init[which_beta]
  }else{
    dat$beta$beta_m <- rep(0, nb)
  }
  if (dat$is_factors && !is.null(init_beta_X_Y)) {
    i_xy <- which(which_beta[, 1] == 1 & which_beta[, 2] == 2)
    dat$beta$beta_m[i_xy] <- init_beta_X_Y
    dat$beta$fix_beta[i_xy] <- TRUE
  }
  dat$beta$beta_s <- rep(0, nb)
  dat$beta$V <- matrix(0, nrow = nb, ncol = nb)

  dat$beta$kl <- 0
  return(dat)
}


