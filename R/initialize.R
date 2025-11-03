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



init_beta <- function(dat, restrict_dag = TRUE, beta_prior_cov = NULL){

  dat$beta <- list()

  B <- check_B_template(dat$B_template, dat$p, restrict_dag = restrict_dag)

  which_beta <- rbind(B$which_tot_u, B$which_tot_c)[,c(2,1), drop=FALSE] ## transpose
  colnames(which_beta) <- c("row", "col")
  dat$beta$beta_j <- which_beta[,1]
  dat$beta$beta_k <- which_beta[,2]
  dat$beta$fix_beta <- c(rep(FALSE, nrow(B$which_tot_u)), rep(TRUE, nrow(B$which_tot_c)))

  nb <- length(dat$beta$beta_j)
  dat$beta$beta_m <- dat$B_init[which_beta]
  dat$beta$beta_s <- rep(0, nb)
  dat$beta$V <- matrix(0, nrow = nb, ncol = nb)

  if (!is.null(beta_prior_cov)) {
    nb <- sum(! dat$beta$fix_beta)
    dat$beta$prior_cov <- check_beta_prior_cov(beta_prior_cov, nb)
    dat$beta$prior_precision <- solve(beta_prior_cov)
  }
  dat$beta$kl <- 0
  return(dat)
}

## Fmatrix will be dat$p by dat$k + 2
## want to estimate first row, except for (1,1) element which is 1
## and second row except for (2,1) which is 0 and (2,2) which is 1
init_beta_factors <- function(dat){
  which_beta <- cbind(c(rep(1,  dat$k - 1),
                        rep(2, dat$k-2)), c(2:dat$k , 3:dat$k))
  colnames(which_beta) <- c("row", "col")
  dat$beta$beta_j <- which_beta[,1]
  dat$beta$beta_k <- which_beta[,2]
  dat$beta$fix_beta <- rep(FALSE, nrow(which_beta))

  nb <- length(dat$beta$beta_j)
  dat$beta$beta_m <- rep(0, nb)
  dat$beta$beta_s <- rep(0, nb)
  dat$beta$V <- matrix(0, nrow = nb, ncol = nb)
  dat$beta$kl <- 0
  return(dat)
}
