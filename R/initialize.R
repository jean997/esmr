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



init_beta <- function(dat, restrict_DAG = TRUE){

  dat$beta <- list()

  B <- check_B_template(dat$B_template, dat$p, restrict_DAG = restrict_DAG)

  which_beta <- rbind(B$which_tot_u, B$which_tot_c)[,c(2,1), drop=FALSE] ## transpose
  colnames(which_beta) <- c("row", "col")
  dat$beta$beta_j <- which_beta[,1]
  dat$beta$beta_k <- which_beta[,2]
  dat$beta$fix_beta <- c(rep(FALSE, nrow(B$which_tot_u)), rep(TRUE, nrow(B$which_tot_c)))

  nb <- length(dat$beta$beta_j)
  dat$beta$beta_m <- dat$B_init[which_beta]
  dat$beta$beta_s <- rep(0, nb)
  dat$beta$V <- matrix(0, nrow = nb, ncol = nb)
  return(dat)
}

