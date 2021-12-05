
#'@title Simulate summary statistics
#'@param Nsnps Number of snps to simulate
#'@param Nx,Ny,Nz sample size for X, Y, and Z_1, ... Z_M gwas. Assuming equal sample size
#'for all possible covariates for simplicity.
#'@param h2x,h2y,h2z Expected heritability of X, Y, and Z_1, ..., Z_M.
#'@param num_conf_z Number of confounders
#'@param num_null_z Number of covariates that don't effect either X or Y
#'@param nsnpx,nsnpy,nsnpz Number of effect SNPs for X, Y, and Zs
#'@param tau tau = h1*(gamma^2)/h2 Proportion of Y heritability explained by X.
#'@param tauz Proportion of Y heritability explained by each confounder. Can be vector of length num_conf_z or scalar.
#'@param omegaz Proportion of X heritability explained by each confounder. Can be vector of length num_conf_z or scalar.
#'@export
sum_stats <- function(nsnps, Nx, Ny, Nz,
                      h2x, h2y, h2z,
                      num_conf_z, num_null_z,
                      nsnpx, nsnpy, nsnpz,
                      taux, tauz, omegaz){



  if(length(tauz) == 1) tauz <- rep(tauz, num_conf_z)
  if(length(omegaz) == 1) omegaz <- rep(omegaz, num_conf_z)

  stopifnot(sum(abs(omegaz)) <1 )
  stopifnot(sum(abs(tauz)) + taux < 1)
  stopifnot(taux >= 0)
  num_z <- num_conf_z + num_null_z
  ## We will assume all SNPs have been normalized to have variance 1 and
  ## draw normalized effects (beta*sqrt(2*f*(1-f))) from a point normal

  #First effect function for Z's
  sigma_z <- sqrt(h2z/nsnpz)
  pz <- nsnpz/nsnps
  gz <- normalmix(pi=c(1-pz, pz),
                  mean=c(0, 0),
                  sd=c(0, sigma_z))

  #Effect function for X
  sigma_x <-sqrt( (1-sum(omegaz))*h2x/nsnpx)
  px <- nsnpx/nsnps
  gx <- normalmix(pi=c(1-px, px),
                  mean=c(0, 0),
                  sd=c(0, sigma_x))


  #Effect function for Y
  sigma_y <-sqrt( (1-(sum(tauz) + taux))*h2y/nsnpy)
  py <- nsnpy/nsnps
  gy <- normalmix(pi=c(1-py, py),
                  mean=c(0, 0),
                  sd=c(0, sigma_y))

  bz <- replicate(n = num_z, rnormalmix(nsnps, gz))
  conf_ix <- seq(num_conf_z)
  null_ix <- seq(num_null_z) + num_conf_z


  gamma_xy <- sqrt(taux*h2y/h2x)
  gamma_zx <- sign(omegaz)*sqrt(abs(omegaz)*h2x/h2z)
  gamma_zy <- sign(tauz)*sqrt(abs(tauz)*h2y/h2z)

  bx <- rnormalmix(nsnps, gx) + bz[, conf_ix] %*% gamma_zx

  by <- rnormalmix(nsnps, gy) + bz[, conf_ix] %*% gamma_zy + gamma_xy*bx


  ### Convert everything from normalized to non-normalized effects
  f <- rbeta(nsnps, 1, 5)
  sebz <- sqrt(1/(2*Nz*f*(1-f)))
  sebx <- sqrt(1/(2*Nx*f*(1-f)))
  seby <- sqrt(1/(2*Ny*f*(1-f)))
  bz <- sebz*sqrt(Nz)*bz
  bx <- bx*sebx*sqrt(Nx)
  by <- by*seby*sqrt(Ny)

  beta_hat_z <- apply(bz, 2, function(b){rnorm(n = nsnps, mean = b, sd = sebz)})
  beta_hat_x <- rnorm(n=nsnps, mean = bx, sd = sebx)
  beta_hat_x2 <- rnorm(n=nsnps, mean = bx, sd = sebx)
  beta_hat_y <- rnorm(n = nsnps, mean = by, sd = seby)

  dat <- cbind(beta_hat_x, beta_hat_x2, sebx, beta_hat_y, seby, beta_hat_z, sebz)
  dat <- data.frame(dat)
  names(dat) <- c("beta_hat_x", "beta_hat_x2", "sebx", "beta_hat_y", "seby",
                  paste0("beta_hat_z", seq(num_z)),
                  "sebz")
  truth <- cbind(bx, by, bz)
  truth <- data.frame(truth)
  names(truth) <- c("bx", "by", paste0("bz", seq(num_z)))
  ret <- list(data = dat, truth = truth, gamma_xy = gamma_xy, gamma_zx = gamma_zx, gamma_zy = gamma_zy)
  return(ret)
}

