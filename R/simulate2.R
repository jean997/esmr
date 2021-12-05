
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
  F_t <- diag(p + 2)

  #confounder
  paYpaX <- which(dir_xz == 1 & dir_yz == 1)
  #mediator
  paYchX <- which(dir_xz == -1 & dir_yz == 1)
  #collider
  chYchX <- which(dir_xz == -1 & dir_yz == -1)

  #parents of x
  paX <- which(dir_xz == 1)
  tot_x_var <- sum(abs(tau_xz[paX]))
  stopifnot(tot_x_var < 1)

  alpha_xz <- rep(0, p)
  alpha_xz[paX] <- sqrt(abs(tau_xz[paX])/(1-tot_x_var))*sign(tau_xz[paX])

  alpha_yz <- rep(0, p)

  #parents of y
  paY <- which(dir_yz == 1)
  tot_y_var <- sum(abs(tau_yz[paY])) + abs(gamma)
  stopifnot(tot_y_var < 1)

  alpha_yz[paX] <- sqrt(abs(tau_yz[paY])/(1-tot_y_var))*sign(tau_yz[paY])
  g <- sqrt(abs(gamma)/(1-tot_y_var))*sign(gamma)

  #parents of each z
  chX <- which(dir_xz == -1)
  tot_z_var <- rep(0, p)
  tot_z_var[paYchX] <- abs(tau_xz[paYchX])
  tot_z_var[chYchX] <- abs(tau_xz[chYchX]) + abs(tau_yz[chYchX])
  alpha_xz[chX] <- sqrt(abs(tau_xz[chX])/(1-tot_z_var[chX]))*sign(tau_xz[chX])
  alpha_yz[chX] <- sqrt(abs(tau_yz[chX])/(1-tot_z_var[chX]))*sign(tau_yz[chX])


  #X to Y
  F_t[1,2] <- g + sum((alpha_xz*alpha_yz)[paYchX])

  #into X
  F_t[paYpaX+2, 1] <- alpha_xz[paYpaX]
  # X to z
  F_t[1,paYchX + 2] <- alpha_xz[paYchX]
  F_t[1, chYchX + 2] <- alpha_xz[chYchX] + g*(alpha_yz[chYchX])

  #Z to Y
  F_t[paYpaX+2, 2] <- (alpha_yz + g*alpha_xz)[paYpaX]
  F_t[paYchX+2,2] <- (alpha_yz)[paYchX]

  #Y to Z
  F_t[2,chYchX + 2] <- alpha_yz[chYchX]

  F_mat <- t(F_t)

  dat <- sim_sumstats_lf(F_mat = F_mat,
                         N = N, J = J, h_2_trait = c(h2_x, h2_y, h2_z),
                         omega = rep(1, p + 2),
                         pi_L = c(pi_x, pi_y, pi_z),
                         overlap_prop = 0,
                         h_2_factor = rep(1, p+2),
                         pi_theta = 1)
  return(dat)
}
