#'@title Simulate multivariate GWAS data
#'@param N Either a single number representing GWAS sample size for all studies or an K+2
#'vector giving sample size for X, Y, and then each of Z_1 to Z_K
#'@param J Number of variants to simulate
#'@param taux_xz,tau_yz Length K vectors. Effects between each of Z_1, ..., Z_K and X and Y respectively.
#'@param dir_xz,dir_yz Direction of effects between each of Z_1, ..., Z_K and X and Y respectively.
#' A value of 1 indicates an effect from Z_k to X. A value of -1 indicates an effect from X to Z_k.
#'@param h2_x,h2_y Scalars giving hertiabiltiy of each of X and Y.
#'@param h2_z Length K vector giving heritability of each of Z_1 to Z_K.
#'@param pi_x,pi_y Scalrs between 0 and 1 giving the proportion of non-zero effects for each of X and Y
#'@param pi_z Length K vector giving proporiton of non-zero effects for each of Z_1 to Z_K.
#'@param gamma Causal effect of X on Y.
#'@details This function simulates data for traits X, Y, and Z_1, ..., Z_k.
#'@export
sim_mv <- function( tau_xz, tau_yz, dir_xz, dir_yz, gamma,
                    h2_x, h2_y, h2_z, N, J,
                    pi_x, pi_y, pi_z){
  p <- length(tau_xz)
  stopifnot(length(tau_yz) == p)
  stopifnot(length(dir_xz) == p)
  stopifnot(length(dir_yz) == p)
  stopifnot(length(h2_z) == p)
  stopifnot(length(pi_z) == p)

  if(any(dir_xz == 1 & dir_yz == -1)){
    stop("No cycles allowed")
  }
  F_t <- matrix(0, nrow = p+2, ncol = p+2)

  #confounder
  paYpaX <- which(dir_xz == 1 & dir_yz == 1)
  #mediator
  paYchX <- which(dir_xz == -1 & dir_yz == 1)
  #collider
  chYchX <- which(dir_xz == -1 & dir_yz == -1)

  #X to Y
  F_t[1,2] <- gamma + sum((tau_xz*tau_yz)[paYchX])

  #into X
  F_t[paYpaX+2, 1] <- tau_xz[paYpaX]
  # X to z
  F_t[1,paYchX + 2] <- tau_xz[paYchX]
  F_t[1, chYchX + 2] <- tau_xz[chYchX] + gamma*(tau_yz[chYchX])

  #Z to Y
  F_t[paYpaX+2, 2] <- (tau_yz + gamma*tau_xz)[paYpaX]
  F_t[paYchX+2,2] <- (tau_yz)[paYchX]

  #Y to Z
  F_t[2,chYchX + 2] <- tau_yz[chYchX]
  d <- 1-colSums(abs(F_t))
  if(any(d <= 0)) stop("Supplied variances are impossible.")
  F_mat <- t(F_t)/d
  F_mat <- sqrt(abs(F_mat))*sign(F_mat)
  diag(F_mat) <- 1

  dat <- sim_sumstats_lf(F_mat = F_mat,
                         N = N, J = J, h_2_trait = c(h2_x, h2_y, h2_z),
                         omega = rep(1, p + 2),
                         pi_L = c(pi_x, pi_y, pi_z),
                         overlap_prop = 0,
                         h_2_factor = rep(1, p+2),
                         pi_theta = 1)
  dat$true_effects <- t(dat$F_mat)/diag(dat$F_mat)
  return(dat)
}
