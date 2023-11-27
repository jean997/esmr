library(GWASBrewer)
library(esmr)

args <- commandArgs(trailingOnly=TRUE)

G <- matrix(
  c(0, 0, 0,
    sqrt(0.4), 0, 0,
    sqrt(0.2), 0, 0),
  nrow = 3,
  byrow = TRUE
)

B_true <- G
B_true[!B_true == 0] <- 1

# Need to trim the cycles
B_alt1 <- matrix(
  c(0, 0, 0,
    1, 0, 0,
    1, 1, 0),
  nrow = 3,
  byrow = TRUE
)

three_sim <- function(seed = NULL) {
  res <- list()
    res$seed <- seed
  true_model <- list()
  full_model <- list()

  h2 <- c(0.5, 0.3, 0.25)
  ## simulate summary statistics
  data(ld_mat_list)
  data(AF)
  set.seed(seed)
  dat <- sim_mv(
    G = G,
    N = 40000,
    J = 5e5,
    h2 = h2,
    pi = 500/5e5,
    R_LD = ld_mat_list,
    af = AF,
    est_s = TRUE
  )

  Z <- with(dat, beta_hat/s_estimate);
  dat$pval <- 2*pnorm(-abs(Z));
  minp <- apply(dat$pval, 1, min)

  # ld pruning
  dat$ld_list_minp <- sim_ld_prune(dat, R_LD = ld_mat_list, pvalue = minp)

  minp <- apply(dat$pval, 1, min)

  ix <- dat$ld_list_minp
  minp <- apply(dat$pval[ix,], 1, min)
  ix1 <- which(minp < 5e-8)

  new_dat <- dat

  true_model <- with(new_dat,
                          esmr(
                            beta_hat_X = beta_hat[ix,],
                            se_X = se_beta_hat[ix,],
                            ix1 = ix1,
                            G = diag(3), # required for network problem
                            direct_effect_template = B_true,
                            max_iter = 300))

  full_model <- with(new_dat,
                          esmr(
                            beta_hat_X = beta_hat[ix,],
                            se_X = se_beta_hat[ix,],
                            ix1 = ix1,
                            G = diag(3), # required for network problem
                            direct_effect_template = B_alt1,
                            max_iter = 300
                          )
  )

  res$true_ll <- logLik(true_model)
  res$true_ll <- logLik(full_model)

  res$log_lik <- lapply(1:3, function(i) {
    pchisq(
      - 2 * (logLik(true_model) - logLik(full_model)),
      df = sum(B_alt1) - sum(B_true),
      lower.tail = FALSE
    )
  })

  res$direct_true <- lapply(true_model, function(x) {
    total_to_direct(t(x$f$fbar) - diag(3))
  })

  res$direct_full <- lapply(full_model, function(x) {
    total_to_direct(t(x$f$fbar) - diag(3))
  })

  return(res)
}

if (length(args) == 0) {
  out_path <- '/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node/'
}
if (length(args) <= 1) {
  seed <- as.integer(Sys.time()) - 1.7e+09
} else {
  seed <- args[[2]]
  out_path <- args[1]
}

saveRDS(
  three_sim(seed = seed),
  file = sprintf('three_node_%s.rds', seed)
  )
