library(simGWAS)
library(esmr)
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

  for(i in 1:3) {
    true_model[[i]] <- with(new_dat,
                            esmr(
                              beta_hat_X = beta_hat[ix,],
                              se_X = se_beta_hat[ix,],
                              ix1 = ix1,
                              G = diag(3), # required for network problem
                              direct_effect_template = B_true,
                              max_iter = 100^i))

    full_model[[i]] <- with(new_dat,
                            esmr(
                              beta_hat_X = beta_hat[ix,],
                              se_X = se_beta_hat[ix,],
                              ix1 = ix1,
                              G = diag(3), # required for network problem
                              direct_effect_template = B_alt1,
                              max_iter = 100^i
                            )
    )

    if (i == 3) break;
    resample_size <- 40000 * 10^i
    cat('resample_size: ', resample_size, '\n')
    new_dat <- resample_sumstats(
      dat = dat, N = resample_size, R_LD = ld_mat_list, af = AF,
      est_s = TRUE
    )

    Z <- with(new_dat, beta_hat/s_estimate);
    new_dat$pval <- 2*pnorm(-abs(Z));

    minp <- apply(new_dat$pval, 1, min)

    # ld pruning
    new_dat$ld_list_minp <- sim_ld_prune(
      new_dat, R_LD = ld_mat_list, pvalue = minp
    )

    minp <- apply(new_dat$pval, 1, min)

    ix <- new_dat$ld_list_minp
    minp <- apply(new_dat$pval[ix,], 1, min)
    ix1 <- which(minp < 5e-8)
  }

  res$log_lik <- lapply(1:3, function(i) {
    pchisq(
      - 2 * (logLik(true_model[[i]]) - logLik(full_model[[i]])),
      df = sum(B_alt1) - sum(B_true),
      lower.tail = FALSE
    )
  })

  res$fbar_true <- lapply(true_model, function(x) {
    total_to_direct(t(x$f$fbar) - diag(3))
  })

  res$fbar_full <- lapply(full_model, function(x) {
    total_to_direct(t(x$f$fbar) - diag(3))
  })
  return(res)
}

sim_results <- list()
# With seed = 6 we get error:
# R_obs is incompatible with trait relationships and heritability.
for (i in 1:5) {
  sim_results[[i]] <- three_sim(seed = i)
}

# sim_results[[6]] <- three_sim(seed = 6)

#length(sim_results)
print(setNames(
  do.call('rbind.data.frame', lapply(sim_results, function(x) x$log_lik)),
  c('40000', '400000', '4000000')))


