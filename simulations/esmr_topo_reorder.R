library(esmr)
library(GWASBrewer)

G <- matrix(c(0, 0, 0, 0, 0,
              sqrt(0.3), 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, -1*sqrt(0.1), 0, 0, 0,
              0, -1*sqrt(0.1), sqrt(0.2), sqrt(0.25), 0), nrow = 5, byrow = 5)

# heritability of each trait
h2 <- c(0.5, 0.3, 0.25, 0.4, 0.3)
## simulate summary statistics
data(ld_mat_list)
data(AF)
set.seed(1)
dat <- sim_mv(G = G,
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

ix <- dat$ld_list_minp
minp <- apply(dat$pval[ix,], 1, min)
ix1 <- which(minp < 5e-8)

# Non lower triangular ordering
inc_order <- c(1,3,5,2,4)

B_template <- dat$direct_trait_effects
B_template[!B_template == 0] <- 1
B_template_reorder <- B_template[inc_order, inc_order]

r2 <- with(dat,
           esmr(beta_hat_X = beta_hat[ix,],
                se_X = se_beta_hat[ix,],
                ix1 = ix1,
                G = diag(5), # required for network problem
                direct_effect_template = B_template))

r2_inc_order <- with(dat,
                     esmr(beta_hat_X = beta_hat[ix,inc_order],
                          se_X = se_beta_hat[ix,inc_order],
                          ix1 = ix1,
                          G = diag(5), # required for network problem
                          direct_effect_template = B_template_reorder))

cat('Direct effects from reordered results:\n')

kableExtra::kable(
  round(esmr:::total_to_direct(t(r2_inc_order$f$fbar) - diag(5)), 3)[r2_inc_order$topo_order, r2_inc_order$topo_order]
  )

stopifnot(
  all(
    t(r2$f$fbar) ==
    t(r2_inc_order$f$fbar)[r2_inc_order$topo_order, r2_inc_order$topo_order]
  )
)


stopifnot(
  all(sort(r2$beta$beta_m) == sort(r2_inc_order$beta$beta_m))
)

stopifnot(
  all(sort(r2$beta$beta_s) == sort(r2_inc_order$beta$beta_s))
)

stopifnot(
  all(sort(r2$beta$beta_s) == sort(r2_inc_order$beta$beta_s))
)
