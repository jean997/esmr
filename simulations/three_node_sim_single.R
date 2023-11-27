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
  res$full_ll <- logLik(full_model)

  res$log_lik <- pchisq(
      - 2 * (res$true_ll - res$full_ll),
      df = sum(B_alt1) - sum(B_true),
      lower.tail = FALSE
    )

  res$true_fbar <- true_model$f$fbar
  res$full_fbar <- full_model$f$fbar

  return(res)
}

n_args <- length(args)
if (n_args == 0) {
  out_path <- '/nfs/turbo/sph-jvmorr/NESMR/simulations/three_node/'
} else {
  out_path <- args[1]
}

if (n_args >= 1) {
  n_seeds <- args[[2]]
}

for (s in 1:n_seeds) {
  seed <- as.integer(Sys.time()) - 1.7e+09
  seed <- as.integer(paste0(sample(1:9, 1), seed))
  cat('Seed set to:', seed, '\nSaving output to: ', out_path, '\n')
  saveRDS(
    three_sim(seed = seed),
    file = file.path(out_path, sprintf('three_node_%s.rds', seed))
  )
}

